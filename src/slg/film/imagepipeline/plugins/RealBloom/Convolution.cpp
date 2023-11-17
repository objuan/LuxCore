#include "Convolution.h"
//#include "ConvolutionThread.h"
#include "ConvolutionFFT.h"

#include <omp.h>
#include <cstring>

void* MemoryManager::memory;
int MemoryManager::maxSize;
int MemoryManager::size;
std::mutex MemoryManager::m_mutex;

void* MemoryManager::_malloc(int m) {
    //  return malloc(m);
    m_mutex.lock();
    if (size + m > maxSize)
    {
   //     SLG_LOG("[RealBloom] RESIZING.. ");

        int newSize = maxSize * 2;
        void* new_memory = malloc(newSize);
        std::memcpy(new_memory, memory, maxSize);
        maxSize = newSize;
        free(memory);
        memory = new_memory;
    }
    void* p = (unsigned char*)memory + size;
    size += m;
    m_mutex.unlock();
    return p;
}

bool _firstTime = true;

namespace RealBloom
{

    void  MyConvolutionThread::start() {
        thread = new std::thread(&MyConvolutionThread::run,this);
    }
    void MyConvolutionThread::run()
    { 
        if (true)
        {
            // Input FFT
           // currStage++;
            //m_status.setFftStage(strFormat("%u/%u %s: Input FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
            fftConv->inputFFT(index);

            //    if (m_status.mustCancel()) throw std::runtime_error();

                // Kernel FFT
         //   currStage++;
          //  m_status.setFftStage(strFormat("%u/%u %s: Kernel FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
            fftConv->kernelFFT(index);

            //  if (m_status.mustCancel()) throw std::runtime_error();

              // Define the name of the arithmetic operation based on deconvolve
         /*   std::string arithmeticName =
                m_capturedParams.methodInfo.FFT_CPU_deconvolve
                ? "Dividing"
                : "Multiplying";*/

                // Multiply/Divide the Fourier transforms
               // currStage++;
              //  m_status.setFftStage(strFormat("%u/%u %s: %s", currStage, numStages, strFromColorChannelID(i).c_str(), arithmeticName.c_str()));
            fftConv->multiplyOrDivide(index);

            //   if (m_status.mustCancel()) throw std::runtime_error();

               // Inverse FFT
            //currStage++;
          //  m_status.setFftStage(strFormat("%u/%u %s: Inverse FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
            fftConv->inverse(index);
        }
        _isEnd = true;

    }

    // ============================================

    static constexpr uint32_t CONV_PROG_TIMESTEP = 1000;

    void drawRect(CmImage* image, int rx, int ry, int rw, int rh);
    void fillRect(CmImage* image, int rx, int ry, int rw, int rh);

    uint32_t ConvolutionStatus::getNumChunksDone() const
    {
        return m_numChunksDone;
    }

    void ConvolutionStatus::setNumChunksDone(uint32_t numChunksDone)
    {
        m_numChunksDone = numChunksDone;
    }

    const std::string& ConvolutionStatus::getFftStage() const
    {
        return m_fftStage;
    }

    void ConvolutionStatus::setFftStage(const std::string& fftStage)
    {
        m_fftStage = fftStage;
    }

    void ConvolutionStatus::reset()
    {
        //super::reset();
        m_numChunksDone = 0;
        m_fftStage = "";
    }

    Convolution::Convolution()
    {}

    ConvolutionParams* Convolution::getParams()
    {
        return &m_params;
    }

    CmImage* Convolution::getImgInputSrc()
    {
        return &m_imgInputSrc;
    }

    void Convolution::setImgInput(CmImage* image)
    {
        m_imgInput = image;
    }

    CmImage* Convolution::getImgKernelSrc()
    {
        return &m_imgKernelSrc;
    }

    void Convolution::setImgKernel(CmImage* image)
    {
        m_imgKernel = image;
    }

    void Convolution::setImgConvPreview(CmImage* image)
    {
        m_imgConvPreview = image;
    }

    void Convolution::setImgConvResult(CmImage* image)
    {
        m_imgConvResult = image;
    }

    void Convolution::previewThreshold(size_t* outNumPixels)
    {
        float threshold = m_params.threshold;
        float transKnee = transformKnee(m_params.knee);
  //      std::atomic_uint64_t numPixels = 0;
        uint64_t numPixels = 0;

        {
            // Input image
         //   std::scoped_lock lock1(*m_imgInput);
            float* inputBuffer = m_imgInput->getImageData();
            uint32_t inputWidth = m_imgInput->getWidth();
            uint32_t inputHeight = m_imgInput->getHeight();

            // Convolution Preview Image
           // std::scoped_lock lock2(*m_imgConvPreview);
            m_imgConvPreview->resize(inputWidth, inputHeight, false);
            float* prevBuffer = m_imgConvPreview->getImageData();

//#pragma omp parallel for
            for (int y = 0; y < inputHeight; y++)
            {
                for (int x = 0; x < inputWidth; x++)
                {
                    uint32_t redIndex = (y * inputWidth + x) * 4;
                    prevBuffer[redIndex + 0] = 0;
                    prevBuffer[redIndex + 1] = 0;
                    prevBuffer[redIndex + 2] = 0;
                    prevBuffer[redIndex + 3] = 1;

                    float v = rgbaToGrayscale(&inputBuffer[redIndex], CONV_THRESHOLD_GRAYSCALE_TYPE);
                    if (v > threshold)
                    {
                        numPixels++;
                        float mul = softThreshold(v, threshold, transKnee);
                        blendAddRGB(prevBuffer, redIndex, inputBuffer, redIndex, mul);
                    }
                }
            }
        }
        m_imgConvPreview->moveToGPU();

        if (outNumPixels)
            *outNumPixels = numPixels;
    }

    void Convolution::previewInput(bool previewMode, std::vector<float>* outBuffer, uint32_t* outWidth, uint32_t* outHeight)
    {
        processInputImage(previewMode, m_params.inputTransformParams, m_imgInputSrc, *m_imgInput, outBuffer, outWidth, outHeight);
    }

    void Convolution::previewKernel(bool previewMode, std::vector<float>* outBuffer, uint32_t* outWidth, uint32_t* outHeight)
    {
        // Return the output dimensions if requested
        if (!previewMode && !outBuffer)
        {
            uint32_t transWidth, transHeight;
            ImageTransform::getOutputDimensions(m_params.kernelTransformParams, m_imgKernelSrc.getWidth(), m_imgKernelSrc.getHeight(), transWidth, transHeight);
            *outWidth = transWidth;
            *outHeight = transHeight;
            return;
        }

        processInputImage(previewMode, m_params.kernelTransformParams, m_imgKernelSrc, *m_imgKernel, outBuffer, outWidth, outHeight);

        bool outerRequest = (!previewMode && outBuffer && outWidth && outHeight);

        // Auto-adjust the exposure
        if (m_params.autoExposure && outerRequest)
        {
            // Get the sum of the grayscale values
            float sumV = 0.0f;
            for (uint32_t y = 0; y < *outHeight; y++)
            {
                for (uint32_t x = 0; x < *outWidth; x++)
                {
                    uint32_t redIndex = (y * *outWidth + x) * 4;

                    float grayscale = rgbaToGrayscale(&(*outBuffer)[redIndex], GrayscaleType::Average);
                    sumV += grayscale;
                }
            }

            // Divide by the sum, and cancel out the convolution multiplier
            if (sumV != 0.0f)
            {
                float mul = 1.0 / ((double)sumV * (double)CONV_MULTIPLIER);
                for (uint32_t i = 0; i < (*outBuffer).size(); i++)
                {
                    if (i % 4 != 3)
                        (*outBuffer)[i] *= mul;
                }
            }
        }
    }

    void Convolution::blend()
    {
        float blendInput = fmaxf(m_params.blendInput, 0.0f);
        float blendConv = fmaxf(m_params.blendConv, 0.0f);
        float blendMix = 0.0f;
        if (!m_params.blendAdditive)
        {
            blendMix = fminf(fmaxf(m_params.blendMix, 0.0f), 1.0f);
            blendConv = blendMix;
            blendInput = 1.0f - blendMix;
        }

        float blendExposure = m_params.blendExposure;

        {
            // Input buffer
          //  std::scoped_lock lock1(m_imgInputCaptured);
            float* inputBuffer = m_imgInputCaptured.getImageData();
            uint32_t inputWidth = m_imgInputCaptured.getWidth();
            uint32_t inputHeight = m_imgInputCaptured.getHeight();

            // Conv. Buffer
           // std::scoped_lock lock2(m_imgOutput);
            float* convBuffer = m_imgOutput.getImageData();
            uint32_t convWidth = m_imgOutput.getWidth();
            uint32_t convHeight = m_imgOutput.getHeight();

            if ((inputWidth == convWidth) && (inputHeight == convHeight))
            {
                // Conv. Result Buffer
              //  std::scoped_lock lock3(*m_imgConvResult);
                m_imgConvResult->resize(inputWidth, inputHeight, false);
                float* convResultBuffer = m_imgConvResult->getImageData();

                float mul = blendConv * getExposureMul(blendExposure);
                for (uint32_t y = 0; y < inputHeight; y++)
                {
                    for (uint32_t x = 0; x < inputWidth; x++)
                    {
                        uint32_t redIndex = (y * inputWidth + x) * 4;

                        convResultBuffer[redIndex + 0] =
                            (inputBuffer[redIndex + 0] * blendInput) + (convBuffer[redIndex + 0] * mul);
                        convResultBuffer[redIndex + 1] =
                            (inputBuffer[redIndex + 1] * blendInput) + (convBuffer[redIndex + 1] * mul);
                        convResultBuffer[redIndex + 2] =
                            (inputBuffer[redIndex + 2] * blendInput) + (convBuffer[redIndex + 2] * mul);

                        convResultBuffer[redIndex + 3] = 1;
                    }
                }
            }
        }
        m_imgConvResult->moveToGPU();
    }

    void Convolution::convolve()
    {
        // If there's already a convolution process going on
        cancel();

        // Capture the parameters to avoid changes during the process
        m_capturedParams = m_params;

        // Reset the status
        m_status.reset();
        //   m_status.setWorking();
        m_status.setFftStage("Initializing");

        // Start the main thread
      /*  m_thread = std::make_shared<std::jthread>([this]()
            {*/
            // Input buffer
        std::vector<float> inputBuffer;
        uint32_t inputWidth = 0, inputHeight = 0;
        previewInput(false, &inputBuffer, &inputWidth, &inputHeight);
        uint32_t inputBufferSize = inputWidth * inputHeight * 4;

        // Kernel buffer
        std::vector<float> kernelBuffer;
        uint32_t kernelWidth = 0, kernelHeight = 0;
        previewKernel(false, &kernelBuffer, &kernelWidth, &kernelHeight);
        uint32_t kernelBufferSize = kernelWidth * kernelHeight * 4;

        convFftCPU(
            kernelBuffer, kernelWidth, kernelHeight,
            inputBuffer, inputWidth, inputHeight, inputBufferSize);


        // std::scoped_lock lock(m_imgInputCaptured);
        m_imgInputCaptured.resize(inputWidth, inputHeight, false);
        float* buffer = m_imgInputCaptured.getImageData();
        std::copy(inputBuffer.data(), inputBuffer.data() + inputBufferSize, buffer);


        // Update convBlendParamsChanged
        /*if (m_status.isOK() && !m_status.mustCancel())
            Async::emitSignal("convBlendParamsChanged", nullptr);

        m_status.setDone();*/
        //  }
     // );
    }

    void Convolution::cancel()
    {
        //if (!m_status.isWorking())
        //{
        //    m_status.reset();
        //    return;
        //}

        //m_status.setMustCancel();

        //if (m_capturedParams.methodInfo.method == ConvolutionMethod::NAIVE_CPU)
        //{
        //    // Tell the sub-threads to stop
        //    for (auto& ct : m_cpuThreads)
        //        ct->stop();
        //}

        //// Wait for the main thread
        //threadJoin(m_thread.get());
        //m_thread = nullptr;

        //m_status.reset();
    }

    const ConvolutionStatus& Convolution::getStatus() const
    {
        return m_status;
    }

    void Convolution::getStatusText(std::string& outStatus, std::string& outMessage, uint32_t& outMessageType) const
    {
        outStatus = "";
        outMessage = "";
        outMessageType = 0;

        //if (m_status.mustCancel())
        //    return;

        //if (!m_status.isOK())
        //{
        //    outMessage = m_status.getError();
        //    outMessageType = 3;
        //}
        //else if (m_status.isWorking())
        //{
        //    float elapsedSec = m_status.getElapsedSec();
        //    if (m_capturedParams.methodInfo.method == ConvolutionMethod::NAIVE_CPU)
        //    {
        //        uint32_t numThreads = m_cpuThreads.size();
        //        uint32_t numPixels = 0;
        //        uint32_t numDone = 0;

        //        for (uint32_t i = 0; i < numThreads; i++)
        //        {
        //            ConvolutionThreadStats* stats = m_cpuThreads[i]->getStats();
        //            if (stats->state == ConvolutionThreadState::Working || stats->state == ConvolutionThreadState::Done)
        //            {
        //                numPixels += stats->numPixels;
        //                numDone += stats->numDone;
        //            }
        //        }

        //        float progress = (numPixels > 0) ? (float)numDone / (float)numPixels : 1.0f;
        //        float remainingSec = (elapsedSec * (float)(numPixels - numDone)) / fmaxf((float)(numDone), EPSILON);
        //        outStatus = strFormat(
        //            "%.1f%%%% (%u/%u)\n%s / %s",
        //            progress * 100.0f,
        //            numDone,
        //            numPixels,
        //            strFromElapsed(elapsedSec).c_str(),
        //            strFromElapsed(remainingSec).c_str());
        //    }
        //    else if (m_capturedParams.methodInfo.method == ConvolutionMethod::NAIVE_GPU)
        //    {
        //        outStatus = strFormat(
        //            "%u/%u chunks\n%s",
        //            m_status.getNumChunksDone(),
        //            m_capturedParams.methodInfo.NAIVE_GPU_numChunks,
        //            strFromElapsed(elapsedSec).c_str());
        //    }
        //    else if (m_capturedParams.methodInfo.method == ConvolutionMethod::FFT_CPU)
        //    {
        //        outStatus = strFormat(
        //            "%s\n%s",
        //            m_status.getFftStage().c_str(),
        //            strFromElapsed(elapsedSec).c_str());
        //    }
        //    else if (m_capturedParams.methodInfo.method == ConvolutionMethod::FFT_GPU)
        //    {
        //        outStatus = strFromElapsed(elapsedSec).c_str();
        //    }
        //}
        //else if (m_status.hasTimestamps())
        //{
        //    float elapsedSec = m_status.getElapsedSec();
        //    outStatus = strFormat("Done (%s)", strFromDuration(elapsedSec).c_str());
        //}
    }

    std::string Convolution::getResourceInfo()
    {
        // Kernel Size
        uint32_t kernelWidth = 0, kernelHeight = 0;
        previewKernel(false, nullptr, &kernelWidth, &kernelHeight);

        // Input Size
        uint32_t inputWidth = m_imgInput->getWidth();
        uint32_t inputHeight = m_imgInput->getHeight();

        // Estimate resource usage
        uint64_t numPixels = 0;
        uint64_t numPixelsPerBlock = 0;
        uint64_t ramUsage = 0;
        uint64_t vramUsage = 0;

        // Number of pixels that pass the threshold
        if ((m_params.methodInfo.method == ConvolutionMethod::NAIVE_CPU) || (m_params.methodInfo.method == ConvolutionMethod::NAIVE_GPU))
            previewThreshold(&numPixels);

        // Number of pixels per thread/chunk
        if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_CPU)
            numPixelsPerBlock = numPixels / m_params.methodInfo.NAIVE_CPU_numThreads;
        else if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_GPU)
            numPixelsPerBlock = numPixels / m_params.methodInfo.NAIVE_GPU_numChunks;

        // Input and kernel size in bytes
        uint64_t inputSizeBytes = (uint64_t)inputWidth * (uint64_t)inputHeight * 4 * sizeof(float);
        uint64_t kernelSizeBytes = (uint64_t)kernelWidth * (uint64_t)kernelHeight * 4 * sizeof(float);

        // Calculate memory usage for different methods
        if (m_params.methodInfo.method == ConvolutionMethod::FFT_CPU)
        {
            std::array<float, 2> kernelOrigin = getKernelOrigin(m_params);
            uint32_t paddedWidth, paddedHeight;
            calcFftConvPadding(
                false, false,
                inputWidth, inputHeight,
                kernelWidth, kernelHeight,
                kernelOrigin[0], kernelOrigin[1],
                paddedWidth, paddedHeight);

            ramUsage = inputSizeBytes + kernelSizeBytes + ((uint64_t)paddedWidth * (uint64_t)paddedHeight * 10 * sizeof(float));
        }
        else if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_CPU)
        {
            ramUsage = inputSizeBytes + kernelSizeBytes + (inputSizeBytes * m_params.methodInfo.NAIVE_CPU_numThreads);
        }
        else if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_GPU)
        {
            // input buffer + kernel buffer + final output buffer + output buffer for chunk
            // + (numPoints * 5) + fbData (same size as input buffer)
            ramUsage = (inputSizeBytes * 4) + kernelSizeBytes + (numPixelsPerBlock * 5 * sizeof(float));

            // (numPoints * 5) + kernel buffer + framebuffer (same size as input buffer)
            vramUsage = (numPixelsPerBlock * 5 * sizeof(float)) + kernelSizeBytes + inputSizeBytes;
        }

        // Format output
      /*  if (m_params.methodInfo.method == ConvolutionMethod::FFT_CPU)
        {
            return strFormat(
                "Est. Memory: %s",
                strFromDataSize(ramUsage).c_str());
        }
        else if (m_params.methodInfo.method == ConvolutionMethod::FFT_GPU)
        {
            return "Undetermined";
        }
        else if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_CPU)
        {
            return strFormat(
                "Total Pixels: %s\nPixels/Thread: %s\nEst. Memory: %s",
                strFromBigInteger(numPixels).c_str(),
                strFromBigInteger(numPixelsPerBlock).c_str(),
                strFromDataSize(ramUsage).c_str());
        }
        else if (m_params.methodInfo.method == ConvolutionMethod::NAIVE_GPU)
        {
            return strFormat(
                "Total Pixels: %s\nPixels/Chunk: %s\nEst. Memory: %s\nEst. VRAM: %s",
                strFromBigInteger(numPixels).c_str(),
                strFromBigInteger(numPixelsPerBlock).c_str(),
                strFromDataSize(ramUsage).c_str(),
                strFromDataSize(vramUsage).c_str());
        }*/

        return "";
    }

    std::array<float, 2> Convolution::getKernelOrigin(const ConvolutionParams& params)
    {
        if (params.useKernelTransformOrigin)
            return params.kernelTransformParams.transform.origin;
        else
            return { 0.5f, 0.5f };
    }

    typedef std::thread*  pThread;

    void Convolution::convFftCPU(
        std::vector<float>& kernelBuffer,
        uint32_t kernelWidth,
        uint32_t kernelHeight,
        std::vector<float>& inputBuffer,
        uint32_t inputWidth,
        uint32_t inputHeight,
        uint32_t inputBufferSize)
    {
      
        try
        {
           // if (_firstTime)
           // {
           //     _firstTime = false;
           //     MemoryManager::Init(400000000);
           // }
           // MemoryManager::Clear();
           //// MemoryManager::Init(400000000);

            //// FFT
            ConvolutionFFT fftConv(
                m_capturedParams, inputBuffer.data(), inputWidth, inputHeight,
                kernelBuffer.data(), kernelWidth, kernelHeight
            );

            uint32_t numStages = 14;
            uint32_t currStage = 0;

            // Prepare padded buffers
            currStage++;
            m_status.setFftStage(strFormat("%u/%u Preparing", currStage, numStages));
            fftConv.pad();

            //  if (m_status.mustCancel()) throw std::runtime_error();

           /* MyConvolutionThread* threads[3];
            for (int i = 0; i < 3; i++)
            {
                threads[i] = new MyConvolutionThread(i,&fftConv);
                threads[i]->start();
            }

            while (true)
            {
                int end = 0;
                for (int i = 0; i < 3; i++)
                    if (threads[i]->isEnd()) end++;
                if (end == 3)
                    break;
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }*/

            // Repeat for 3 color channels
            for (uint32_t i = 0; i < 3; i++)
            {
                // Input FFT
                currStage++;
                m_status.setFftStage(strFormat("%u/%u %s: Input FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
                fftConv.inputFFT(i);

                //    if (m_status.mustCancel()) throw std::runtime_error();

                    // Kernel FFT
                currStage++;
                m_status.setFftStage(strFormat("%u/%u %s: Kernel FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
                fftConv.kernelFFT(i);

                //  if (m_status.mustCancel()) throw std::runtime_error();

                  // Define the name of the arithmetic operation based on deconvolve
                std::string arithmeticName =
                    m_capturedParams.methodInfo.FFT_CPU_deconvolve
                    ? "Dividing"
                    : "Multiplying";

                // Multiply/Divide the Fourier transforms
                currStage++;
                m_status.setFftStage(strFormat("%u/%u %s: %s", currStage, numStages, strFromColorChannelID(i).c_str(), arithmeticName.c_str()));
                fftConv.multiplyOrDivide(i);

                //   if (m_status.mustCancel()) throw std::runtime_error();

                   // Inverse FFT
                currStage++;
                m_status.setFftStage(strFormat("%u/%u %s: Inverse FFT", currStage, numStages, strFromColorChannelID(i).c_str()));
                fftConv.inverse(i);

                // if (m_status.mustCancel()) throw std::runtime_error();
            }

            //for (int i = 0; i < 3; i++)
            //     delete threads[i];
         
            // Get the final output
            currStage++;
            m_status.setFftStage(strFormat("%u/%u Finalizing", currStage, numStages));
            fftConv.output();

            // Update the output image
            {
                // std::scoped_lock lock(m_imgOutput);
                m_imgOutput.resize(inputWidth, inputHeight, false);
                float* convBuffer = m_imgOutput.getImageData();
                std::copy(
                    fftConv.getBuffer().data(),
                    fftConv.getBuffer().data() + fftConv.getBuffer().size(),
                    convBuffer
                );
            }
            //MemoryManager::Dispose();
        }
        catch (const std::runtime_error& e)
        {
            //  m_status.setError(e.what());
        }
    }


    void drawRect(CmImage* image, int rx, int ry, int rw, int rh)
    {
      //  std::scoped_lock lock(*image);
        float* buffer = image->getImageData();
        uint32_t imageWidth = image->getWidth();
        uint32_t imageHeight = image->getHeight();

        int rx2 = rx + rw;
        int ry2 = ry + rh;

        uint32_t redIndex;
        for (int y = ry; y < ry2; y++)
        {
            if (checkBounds(rx, y, imageWidth, imageHeight))
            {
                redIndex = (y * imageWidth + rx) * 4;
                buffer[redIndex + 0] = 0;
                buffer[redIndex + 1] = 0;
                buffer[redIndex + 2] = 1;
                buffer[redIndex + 3] = 1;
            }

            if (checkBounds(rx2, y, imageWidth, imageHeight))
            {
                redIndex = (y * imageWidth + rx2) * 4;
                buffer[redIndex + 0] = 0;
                buffer[redIndex + 1] = 0;
                buffer[redIndex + 2] = 1;
                buffer[redIndex + 3] = 1;
            }
        }
        for (int x = rx; x < rx2; x++)
        {
            if (checkBounds(x, ry, imageWidth, imageHeight))
            {
                redIndex = (ry * imageWidth + x) * 4;
                buffer[redIndex + 0] = 0;
                buffer[redIndex + 1] = 0;
                buffer[redIndex + 2] = 1;
                buffer[redIndex + 3] = 1;
            }

            if (checkBounds(x, ry2, imageWidth, imageHeight))
            {
                redIndex = (ry2 * imageWidth + x) * 4;
                buffer[redIndex + 0] = 0;
                buffer[redIndex + 1] = 0;
                buffer[redIndex + 2] = 1;
                buffer[redIndex + 3] = 1;
            }
        }
    }

    void fillRect(CmImage* image, int rx, int ry, int rw, int rh)
    {
       // std::scoped_lock lock(*image);
        float* buffer = image->getImageData();
        uint32_t imageWidth = image->getWidth();
        uint32_t imageHeight = image->getHeight();

        int rx2 = rx + rw;
        int ry2 = ry + rh;

        for (int y = ry; y < ry2; y++)
        {
            for (int x = rx; x < rx2; x++)
            {
                if (checkBounds(x, y, imageWidth, imageHeight))
                {
                    uint32_t redIndex = (y * imageWidth + x) * 4;
                    buffer[redIndex + 0] = 0;
                    buffer[redIndex + 1] = 0;
                    buffer[redIndex + 2] = 1;
                    buffer[redIndex + 3] = 1;
                }
            }
        }
    }

}
