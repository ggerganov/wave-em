/*! \file main.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "build_timestamp.h"

#include "fftw3.h"
#include "reed-solomon/rs.hpp"

#include <SDL2/SDL.h>
#include <SDL2/SDL_audio.h>

#include <cmath>
#include <cstdio>
#include <array>
#include <string>
#include <chrono>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#ifdef __EMSCRIPTEN__
#include "emscripten/emscripten.h"
#endif

#ifdef main
#undef main
#endif

struct DataRxTx;

static SDL_AudioDeviceID devid_in = 0;
static SDL_AudioDeviceID devid_out = 0;

std::shared_ptr<DataRxTx> g_data;

namespace Data {
    namespace Constants {
        constexpr auto kSubFrames = 32;
        constexpr auto kMaxSamplesPerFrame = 1024;
        constexpr auto kMaxDataBits = 256;
        constexpr auto kMaxBitsPerChecksum = 2;
        constexpr auto kMaxSpectrumHistory = 2;
        constexpr auto ikMaxSpectrumHistory = 1.0/kMaxSpectrumHistory;
        constexpr auto kMaxDataSize = 1024;
    }

    using AmplitudeData = std::array<float, 2*Constants::kMaxSamplesPerFrame>;
    using SpectrumData  = std::array<float,   Constants::kMaxSamplesPerFrame>;
}

namespace {
    constexpr float IRAND_MAX = 1.0f/RAND_MAX;
    inline float frand() { return ((float)(rand()%RAND_MAX)*IRAND_MAX); }

    inline void addAmplitude(const ::Data::AmplitudeData & src, ::Data::AmplitudeData & dst, float scalar, int startId, int finalId) {
        for (int i = startId; i < finalId; i++) {
            dst[i] += scalar*src[i];
        }
    }

    template <class T>
        float getTime_ms(const T & tStart, const T & tEnd) {
            return ((float)(std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count()))/1000.0;
        }
}

struct DataRxTx {
    bool needUpdate = false;

    int paramFreqDelta = 3;
    int paramFreqStart = 140;
    int paramFreqCheck = 342;
    int paramFramesPerTx = 4;
    int paramBytesPerTx = 8;
    int paramECCBytesPerTx = 4;
    int paramVolume = 20;

    DataRxTx(int aSampleRate, int aSamplesPerFrame, int aSampleSizeB, const char * text) {
        sampleSizeBytes = aSampleSizeB;
        sampleRate = aSampleRate;
        samplesPerFrame = aSamplesPerFrame;

        init(text);
    }

    void init(const char * text) {
        nIterations = 0;
        encodeIdParity = true;
        hasData = false;

        samplesPerSubFrame = (samplesPerFrame/::Data::Constants::kSubFrames);
        isamplesPerFrame = 1.0f/samplesPerFrame;

        sendVolume = ((float)(paramVolume))/100.0f;
        hzPerFrame = sampleRate/samplesPerFrame;
        ihzPerFrame = 1.0/hzPerFrame;

        freqDelta_hz = paramFreqDelta*hzPerFrame;
        freqStart_hz = paramFreqStart*hzPerFrame;
        freqCheck_hz = paramFreqCheck*hzPerFrame;

        frameId = 0;
        nRampFrames = 64;
        nRampFramesBegin = 64;
        nRampFramesEnd = 64;
        nRampFramesBlend = 64;
        dataId = 0;

        sendId = 0;
        subFramesPerTx = paramFramesPerTx*::Data::Constants::kSubFrames;
        curTxSubFrameId = 0;
        nDataBitsPerTx = paramBytesPerTx*8;
        nECCBytesPerTx = paramECCBytesPerTx;

        outputBlock.fill(0);
        outputBlockTmp.fill(0);
        sendData.fill(0);

        textToSend = std::string(text);

        for (int k = 0; k < (int) dataBits.size(); ++k) {
            auto freq = freqStart_hz + freqDelta_hz*k;
            dataFreqs_hz[k] = freq;

            float phaseOffset = 2*M_PI*::frand();
            for (int i = 0; i < samplesPerFrame; i++) {
                bitAmplitude[k][i] = std::sin((2.0*M_PI*i)*freq*isamplesPerFrame*ihzPerFrame + phaseOffset);
            }
            for (int i = 0; i < samplesPerFrame; i++) {
                bit0Amplitude[k][i] = std::sin((2.0*M_PI*i)*(freq + hzPerFrame)*isamplesPerFrame*ihzPerFrame + phaseOffset);
            }

            //printf("Bit %d -> %4.2f Hz\n", k, freq);
        }

        for (int k = 0; k < ::Data::Constants::kMaxBitsPerChecksum; ++k) {
            checksumAmplitude[k].fill(0);
            checksum0Amplitude[k].fill(0);
        }

        for (int k = 0; k < ::Data::Constants::kMaxBitsPerChecksum; ++k) {
            auto freq = freqCheck_hz + freqDelta_hz*k;
            checksumFreqs_hz[k] = freq;

            float phaseOffset = 2*M_PI*::frand();
            for (int i = 0; i < samplesPerFrame; i++) {
                checksumAmplitude[k][i] = std::sin((2.0*M_PI*i)*freq*isamplesPerFrame*ihzPerFrame + phaseOffset);
            }
            for (int i = 0; i < samplesPerFrame; i++) {
                checksum0Amplitude[k][i] = std::sin((2.0*M_PI*i)*(freq + hzPerFrame)*isamplesPerFrame*ihzPerFrame + phaseOffset);
            }
        }

        if (nECCBytesPerTx >= nDataBitsPerTx/8) {
            rs.reset();
            nECCBytesPerTx = 0;
        } else {
            rs = std::make_shared<RS::ReedSolomon>(nDataBitsPerTx/8 - nECCBytesPerTx, nECCBytesPerTx);
        }

        for (int i = 0; i < (int) textToSend.length(); ++i) {
            sendData[i] = textToSend[i];
        }
        sendData[textToSend.length()] = 0;
        if (textToSend.length() > 0) {
            hasData = true;
        }

        // Rx
        receivingData = false;

        historyId = 0;
        receivedId = 0;
        nConfirmFrames = 16;

        sampleAmplitude.fill(0);

        sampleSpectrum.fill(0);
        sampleSpectrumTmp.fill(0);
        historySpectrumAverage.fill(0);

        for (auto & s : historySpectrum) {
            s.fill(0);
        }
        rxData.fill(0);

        if (fftPlan) fftwf_destroy_plan(fftPlan);
        if (fftIn) fftwf_free(fftIn);
        if (fftOut) fftwf_free(fftOut);

        fftIn = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*samplesPerFrame);
        fftOut = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*samplesPerFrame);
        fftPlan = fftwf_plan_dft_r2c_1d(1*samplesPerFrame, sampleAmplitude.data(), fftOut, FFTW_ESTIMATE);

        for (int i = 0; i < samplesPerFrame; ++i) {
            fftOut[i][0] = 0.0f;
            fftOut[i][1] = 0.0f;
        }

        //printf("Initialized succefully\n");
    }

    void send() {
        frameId = 0;
        nIterations = 0;

        //{
        //    for (int iTrigger = 0; iTrigger < 16; ++iTrigger) {
        //        outputBlock.fill(0);
        //        ::addAmplitude(checksumAmplitude[0], outputBlock, 0.1*sendVolume, 0, samplesPerFrame);
        //        SDL_QueueAudio(devid_out, outputBlock.data(), samplesPerFrame*sampleSizeBytes);
        //    }
        //}

        while(hasData) {
            for (int sfi = 0; sfi < ::Data::Constants::kSubFrames; ++sfi) {
                auto subFrame = nIterations % ::Data::Constants::kSubFrames;

                auto sampleStartId = (subFrame*samplesPerSubFrame);
                auto sampleFinalId = sampleStartId + samplesPerSubFrame;

                if (curTxSubFrameId >= subFramesPerTx) {
                    curTxSubFrameId = 0;
                    frameId = 0;
                    sendId += nDataBitsPerTx/8 - nECCBytesPerTx;
                } else if (curTxSubFrameId >= nRampFrames) {
                    nRampFrames = nRampFramesBlend;
                }

                if (sendData[sendId] == 0) {
                    nRampFrames = nRampFramesEnd;
                } else {
                    curTxSubFrameId = frameId;

                    static std::array<std::uint8_t, ::Data::Constants::kMaxDataBits/8> encoded;
                    if (rs) {
                        rs->Encode(sendData.data() + sendId, encoded.data());
                    } else {
                        for (int j = 0; j < nDataBitsPerTx/8; ++j) {
                            encoded[j] = sendData[sendId + j];
                        }
                    }

                    for (int j = 0; j < nDataBitsPerTx/8; ++j) {
                        for (int i = 0; i < 8; ++i) {
                            dataBits[j*8 + i] = encoded[j] & (1 << i);
                        }
                    }
                }

                std::uint16_t nFreq = 0;
                std::uint16_t checksum = 0;

                checksum += (1 << 0);

                if (encodeIdParity) {
                    if ((dataId + sendId/(nDataBitsPerTx/8 - nECCBytesPerTx)) & 1) {
                        checksum += (1 << 1);
                    }
                }

                for (int i = sampleStartId; i < sampleFinalId; ++i) {
                    outputBlockTmp[i] = 0.0f;
                }

                for (int k = 0; k < nDataBitsPerTx; ++k) {
                    ++nFreq;
                    if (dataBits[k] == false) {
                        checksum += (1 << ((k%8)+2));
                        ::addAmplitude(bit0Amplitude[k], outputBlockTmp, sendVolume, sampleStartId, sampleFinalId);
                        continue;
                    }
                    ::addAmplitude(bitAmplitude[k], outputBlockTmp, sendVolume, sampleStartId, sampleFinalId);
                }

                for (int k = 0; k < ::Data::Constants::kMaxBitsPerChecksum; ++k) {
                    ++nFreq;
                    if ((checksum & (1 << k)) || (k == 0)) {
                        ::addAmplitude(checksumAmplitude[k], outputBlockTmp, sendVolume, sampleStartId, sampleFinalId);
                        continue;
                    }
                    ::addAmplitude(checksum0Amplitude[k], outputBlockTmp, sendVolume, sampleStartId, sampleFinalId);
                }

                if (nFreq == 0) nFreq = 1;
                float scale = 1.0f/nFreq;
                for (int i = sampleStartId; i < sampleFinalId; ++i) {
                    outputBlockTmp[i] *= scale;
                }

                static double interp = 0.0;
                if (frameId == 0 && sendId == 0) {
                    interp = 0.0f;
                }
                double dinterp = 1.0/(nRampFrames*samplesPerSubFrame);
                if (frameId < nRampFrames) {
                    for (int i = sampleStartId; i < sampleFinalId; ++i) {
                        interp = std::min(1.0, interp + dinterp);
                        outputBlock[i] = interp*outputBlockTmp[i];
                    }
                } else if (subFramesPerTx > 0 && frameId >= subFramesPerTx - nRampFrames) {
                    for (int i = sampleStartId; i < sampleFinalId; ++i) {
                        interp = std::max(0.0, interp - dinterp);
                        outputBlock[i] = interp*outputBlockTmp[i];
                    }
                } else {
                    interp = 1.0;
                    for (int i = sampleStartId; i < sampleFinalId; ++i) {
                        outputBlock[i] = outputBlockTmp[i];
                    }
                }


                if (sendId >= textToSend.length() && interp == 0.0) {
                    textToSend = "";
                    hasData = false;
                }

                ++frameId;
                ++nIterations;
            }

            SDL_QueueAudio(devid_out, outputBlock.data(), samplesPerFrame*sampleSizeBytes);
        }

        frameId = 0;
        nIterations = 0;
    }

    void receive() {
        static int nCalls = 0;
        static float tSum_ms = 0.0f;
        auto tCallStart = std::chrono::high_resolution_clock::now();

        if (needUpdate) {
            init("");
            needUpdate = false;
        }

        while (hasData == false) {
            auto subFrame = nIterations % ::Data::Constants::kSubFrames;

            auto sampleStartId = (subFrame*samplesPerSubFrame);
            auto sampleFinalId = sampleStartId + samplesPerSubFrame;

            // check if receiving data
            {
                static std::array<std::uint8_t, ::Data::Constants::kMaxDataBits/8> receivedData;
                static std::array<std::uint8_t, ::Data::Constants::kMaxDataBits/8> receivedDataLast;

                std::uint16_t requiredChecksum = 0;
                std::uint16_t curChecksum = 0;
                std::uint8_t curParity = 0;

                static std::uint8_t lastParity = 2;
                static std::uint16_t lastChecksum = -1;
                static std::uint16_t lastReceivedChecksum = -1;
                static std::uint16_t nTimesReceived = 0;

                if (receivedId == 0) {
                    receivedDataLast.fill(0);
                    lastReceivedChecksum = 0;
                }

                receivedData.fill(0);
                requiredChecksum += 1;

                bool isValid = true;
                {
                    int bin = std::round(checksumFreqs_hz[0]*ihzPerFrame);
                    if (historySpectrumAverage[bin] < 10*historySpectrumAverage[bin - 1] &&
						historySpectrumAverage[bin] < 10*historySpectrumAverage[bin + 1]) {
                        if (receivingData == true) {
                            receivingData = false;
                        }
                    } else {
                        curChecksum += 1;
                        if (receivingData == false) {
                            receivingData = true;
                        }
                    }
                }

                for (int k = 0; k < nDataBitsPerTx; ++k) {
                    int bin = std::round(dataFreqs_hz[k]*ihzPerFrame);
                    if (historySpectrumAverage[bin] > 1.0*historySpectrumAverage[bin + 1]) {
                        receivedData[k/8] += (1 << (k%8));
                    }
                }

                for (int k = 1; k < ::Data::Constants::kMaxBitsPerChecksum; ++k) {
                    int bin = std::round(checksumFreqs_hz[k]*ihzPerFrame);
                    if (historySpectrumAverage[bin] > 1.0*historySpectrumAverage[bin + 1]) {
                        curChecksum += (1 << k);
                        if (k == 1) curParity = 1;
                    }
                }

                requiredChecksum = (requiredChecksum & ((1 << ::Data::Constants::kMaxBitsPerChecksum) - 1));

                isValid = receivingData;
                bool checksumMatch = (lastChecksum == curChecksum);

                if (rs) {
                    static std::array<std::uint8_t, ::Data::Constants::kMaxDataBits/8> repaired;
                    bool decoded = true;
                    if (rs->Decode(receivedData.data(), repaired.data()) != 0) {
                        decoded = false;
                    } else {
                        for (int i = 0; i < nDataBitsPerTx/8 - nECCBytesPerTx; ++i) {
                            receivedData[i] = repaired[i];
                        }
                        receivedData[nDataBitsPerTx/8 - nECCBytesPerTx] = 0;
                        if (receivedData[0] == 0) decoded = false;
                    }
                    checksumMatch = true;
                    isValid &= decoded;
                }

                if (isValid && checksumMatch) {
                    for (int i = 0; i < nDataBitsPerTx/8 - nECCBytesPerTx; ++i) {
                        if (receivedData[i] == 0) receivedData[i] = ' ';
                    }
                    if (++nTimesReceived == nConfirmFrames && receivedData != receivedDataLast) {
                        receivedDataLast = receivedData;
                        lastReceivedChecksum = curChecksum;

                        //printf("Receiving data: %s\n", receivedData.data());
                        static auto tLast = std::chrono::steady_clock::now();
                        auto tNow = std::chrono::steady_clock::now();
                        if (getTime_ms(tLast, tNow) > 500) {
                            receivedId = 0;
                            rxData.fill(0);
                        } else {
                            if (curParity == lastParity && receivedId > 0 && encodeIdParity) {
                                receivedId -= nDataBitsPerTx/8 - nECCBytesPerTx;
                            }
                        }
                        lastParity = curParity;
                        tLast = tNow;

                        for (int i = 0; i < nDataBitsPerTx/8 - nECCBytesPerTx; ++i) {
                            rxData[receivedId++] = receivedData[i];
                        }
                    }
                } else if (isValid && (checksumMatch == false)) {
                    lastChecksum = curChecksum;
                    nTimesReceived = 0;
                } else if (isValid == false) {
                    lastChecksum = -1;
                    nTimesReceived = 0;
                }
            }

            // read capture data

            SDL_DequeueAudio(devid_in, sampleAmplitude.data() + sampleStartId, samplesPerSubFrame*sampleSizeBytes);

            if (sampleStartId == 0) {
                // calculate spectrum
                for (int i = 0; i < samplesPerFrame; ++i) {
                    fftIn[i][0] = sampleAmplitude[i];
                    fftIn[i][1] = 0;
                }

                fftwf_execute(fftPlan);

                for (int i = 0; i < samplesPerFrame; ++i) {
                    sampleSpectrumTmp[i] = (fftOut[i][0]*fftOut[i][0] + fftOut[i][1]*fftOut[i][1]);
                }
                for (int i = 1; i < samplesPerFrame/2; ++i) {
                    sampleSpectrumTmp[i] += sampleSpectrumTmp[samplesPerFrame - i];
                    sampleSpectrumTmp[samplesPerFrame - i] = 0.0f;
                }

                sampleSpectrum = sampleSpectrumTmp;

                // store spectrum in history
                static int nNotReceiving = 0;
                if (receivingData == false) ++nNotReceiving; else nNotReceiving = 0;
                if (nNotReceiving == 8*::Data::Constants::kSubFrames) {
                    for (auto & s : historySpectrum) {
                        s.fill(0);
                    }
                    historySpectrumAverage.fill(0);
                }

                if (true) {
                    for (int i = 0; i < samplesPerFrame/2; ++i) {
                        historySpectrumAverage[i] *= ::Data::Constants::kMaxSpectrumHistory;
                        historySpectrumAverage[i] -= historySpectrum[historyId][i];
                        historySpectrumAverage[i] += sampleSpectrum[i];
                        historySpectrumAverage[i] *= ::Data::Constants::ikMaxSpectrumHistory;
                    }
                    historySpectrum[historyId] = sampleSpectrum;
                    if (++historyId >= ::Data::Constants::kMaxSpectrumHistory) historyId = 0;
                }
            }

            ++nIterations;

            if (nIterations % 1024 == 0) {
                break;
            }
        }

        auto tCallEnd = std::chrono::high_resolution_clock::now();
        tSum_ms += getTime_ms(tCallStart, tCallEnd);
        if (++nCalls == 10) {
            averageRxTime_ms = tSum_ms/nCalls;
            tSum_ms = 0.0f;
            nCalls = 0;
        }

        if (SDL_GetQueuedAudioSize(devid_in) > 0) {
            printf("nIter = %d, Queue size: %d\n", nIterations, SDL_GetQueuedAudioSize(devid_in));
            if (SDL_GetQueuedAudioSize(devid_in) > 128*1024) {
                SDL_ClearQueuedAudio(devid_in);
            }
        }
    }

    int nIterations;
    bool encodeIdParity;

    // Rx
    bool receivingData;

    int historyId;
    int receivedId;
    int nConfirmFrames;

    fftwf_plan fftPlan = 0;
    fftwf_complex *fftIn = 0;
    fftwf_complex *fftOut = 0;

    ::Data::AmplitudeData sampleAmplitude;

    ::Data::SpectrumData sampleSpectrum;
    ::Data::SpectrumData sampleSpectrumTmp;
    ::Data::SpectrumData historySpectrumAverage;

    std::array<::Data::SpectrumData, ::Data::Constants::kMaxSpectrumHistory> historySpectrum;
    std::array<char, ::Data::Constants::kMaxDataSize> rxData;

    // Tx
    bool hasData;
    int sampleSizeBytes;
    int sampleRate;
    int samplesPerFrame;
    int samplesPerSubFrame;
    float isamplesPerFrame;

    ::Data::AmplitudeData outputBlock;
    ::Data::AmplitudeData outputBlockTmp;

    std::array<::Data::AmplitudeData, ::Data::Constants::kMaxDataBits> bitAmplitude;
    std::array<::Data::AmplitudeData, ::Data::Constants::kMaxDataBits> bit0Amplitude;
    std::array<::Data::AmplitudeData, ::Data::Constants::kMaxBitsPerChecksum> checksumAmplitude;
    std::array<::Data::AmplitudeData, ::Data::Constants::kMaxBitsPerChecksum> checksum0Amplitude;

    float sendVolume;
    float hzPerFrame;
    float ihzPerFrame;

    float freqStart_hz;
    float freqDelta_hz;
    float freqCheck_hz;

    int frameId;
    int nRampFrames;
    int nRampFramesBegin;
    int nRampFramesEnd;
    int nRampFramesBlend;
    int dataId;

    std::array<bool, ::Data::Constants::kMaxDataBits> dataBits;
    std::array<float, ::Data::Constants::kMaxDataBits> dataFreqs_hz;
    std::array<float, ::Data::Constants::kMaxBitsPerChecksum> checksumFreqs_hz;

    int sendId;
    int subFramesPerTx;
    int curTxSubFrameId;
    int nDataBitsPerTx;
    int nECCBytesPerTx;

    std::array<char, ::Data::Constants::kMaxDataSize> sendData;

    std::shared_ptr<RS::ReedSolomon> rs = nullptr;

    float averageRxTime_ms = 0.0;

    std::string textToSend;
};

extern "C" {
    int setText(const char * text) {
        g_data->init(text);
        return 0;
    }

    const char * getText() {
        return g_data->rxData.data();
    }

    int getSampleRate() {
        return g_data->sampleRate;
    }

    float getAverageRxTime_ms() {
        return g_data->averageRxTime_ms;
    }

    void setParameters(
        int paramFreqDelta,
        int paramFreqStart,
        int paramFreqCheck,
        int paramFramesPerTx,
        int paramBytesPerTx,
        int paramECCBytesPerTx,
        int paramVolume) {
        if (g_data == nullptr) return;

        //printf("Params: %d %d %d %d %d %d %d\n",
        //       paramFreqDelta,
        //       paramFreqStart,
        //       paramFreqCheck,
        //       paramFramesPerTx,
        //       paramBytesPerTx,
        //       paramECCBytesPerTx,
        //       paramVolume);

        g_data->paramFreqDelta = paramFreqDelta;
        g_data->paramFreqStart = paramFreqStart;
        g_data->paramFreqCheck = paramFreqCheck;
        g_data->paramFramesPerTx = paramFramesPerTx;
        g_data->paramBytesPerTx = paramBytesPerTx;
        g_data->paramECCBytesPerTx = paramECCBytesPerTx;
        g_data->paramVolume = paramVolume;

        g_data->needUpdate = true;
    }
}

void update() {
    SDL_Event e;

    while (SDL_PollEvent(&e)) { }

    if (g_data->hasData == false) {
        SDL_PauseAudioDevice(devid_out, SDL_FALSE);
        SDL_PauseAudioDevice(devid_in, SDL_FALSE);

        static auto tLastNoData = std::chrono::high_resolution_clock::now();
        auto tNow = std::chrono::high_resolution_clock::now();

        if (SDL_GetQueuedAudioSize(devid_out) == 0) {
            if (::getTime_ms(tLastNoData, tNow) > 500.0f) {
                g_data->receive();
            } else {
                SDL_ClearQueuedAudio(devid_in);
            }
        } else {
            tLastNoData = tNow;
            SDL_ClearQueuedAudio(devid_in);
        }
    } else {
        SDL_PauseAudioDevice(devid_out, SDL_TRUE);
        SDL_PauseAudioDevice(devid_in, SDL_TRUE);

        g_data->send();
    }
}

int main(int argc, char** argv) {
    printf("Build time: %s\n", BUILD_TIMESTAMP);

    const char *captureDeviceName = argv[1];

    SDL_LogSetPriority(SDL_LOG_CATEGORY_APPLICATION, SDL_LOG_PRIORITY_INFO);

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return (1);
    }

    {
        int devcount = SDL_GetNumAudioDevices(SDL_FALSE);
        for (int i = 0; i < devcount; i++) {
            printf("Output  device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_FALSE));
        }
    }
    {
        int devcount = SDL_GetNumAudioDevices(SDL_TRUE);
        for (int i = 0; i < devcount; i++) {
            printf("Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
        }
    }

    SDL_AudioSpec desiredSpec;
    SDL_zero(desiredSpec);

    desiredSpec.freq = 44100;
    desiredSpec.format = AUDIO_F32SYS;
    desiredSpec.channels = 1;
    desiredSpec.samples = 1024;
    desiredSpec.callback = NULL;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    devid_out = SDL_OpenAudioDevice(NULL, SDL_FALSE, &desiredSpec, &obtainedSpec, SDL_AUDIO_ALLOW_ANY_CHANGE);
    if (!devid_out) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for playback: %s!\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }

    printf("Obtained spec for output device (SDL Id = %d):\n", devid_out);
    printf("    - Sample rate:       %d (make sure this is the same between devices)\n", obtainedSpec.freq);
    printf("    - Format:            %d (required: %d)\n", obtainedSpec.format, desiredSpec.format);
    printf("    - Channels:          %d (required: %d)\n", obtainedSpec.channels, desiredSpec.channels);
    printf("    - Samples per frame: %d (required: %d)\n", obtainedSpec.samples, desiredSpec.samples);

    if (obtainedSpec.format != desiredSpec.format ||
        obtainedSpec.channels != desiredSpec.channels ||
        obtainedSpec.samples != desiredSpec.samples) {
        SDL_CloseAudio();
        throw std::runtime_error("Failed to initialize desired SDL_OpenAudio!");
    }

    SDL_AudioSpec captureSpec;
    captureSpec = obtainedSpec;

    SDL_Log("Opening capture device %s%s%s...\n",
            captureDeviceName ? "'" : "",
            captureDeviceName ? captureDeviceName : "[[default]]",
            captureDeviceName ? "'" : "");

    devid_in = SDL_OpenAudioDevice(argv[1], SDL_TRUE, &captureSpec, &captureSpec, 0);
    if (!devid_in) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for capture: %s!\n", SDL_GetError());
        SDL_Quit();
        exit(1);
    }

    int sampleSizeBytes = 4;
    switch (obtainedSpec.format) {
        case AUDIO_U8:
        case AUDIO_S8:
            sampleSizeBytes = 1;
            break;
        case AUDIO_U16SYS:
        case AUDIO_S16SYS:
            sampleSizeBytes = 2;
            break;
        case AUDIO_S32SYS:
        case AUDIO_F32SYS:
            sampleSizeBytes = 4;
            break;
    }

    g_data = std::make_shared<DataRxTx>(
        obtainedSpec.freq,
        obtainedSpec.samples,
        sampleSizeBytes,
        "");

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(update, 0, 1);
#else
    while(true) {
        SDL_Delay(20);
        update();
    }
#endif

    return 0;
}
