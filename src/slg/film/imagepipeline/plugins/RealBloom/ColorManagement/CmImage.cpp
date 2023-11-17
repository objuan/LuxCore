#include "CmImage.h"

#include <mutex>
#include "CMS.h"

//
//template <typename T>
//void clearVector(std::vector<T>& v)
//{
//    v.clear();
//    std::vector<T>().swap(v);
//}


//std::shared_ptr<GlFramebuffer> CmImage::s_framebuffer = nullptr;

CmImage::CmImage(const std::string& id, const std::string& name, uint32_t width, uint32_t height, std::array<float, 4> fillColor, bool useExposure, bool useGlobalFB)
    : m_id(id), m_name(name), m_width(width), m_height(height), m_useExposure(useExposure), m_useGlobalFB(useGlobalFB)
{
    lock();
    resize(std::max(width, 1u), std::max(height, 1u), false);
    fill(fillColor, false);
    unlock();
}

CmImage::~CmImage()
{
   // std::scoped_lock lock(m_mutex);
}

const std::string& CmImage::getID() const
{
    return m_id;
}

const std::string& CmImage::getName() const
{
    return m_name;
}

const std::string& CmImage::getSourceName() const
{
    return m_sourceName;
}

void CmImage::setSourceName(const std::string& sourceName)
{
    m_sourceName = sourceName;
}

uint32_t CmImage::getWidth() const
{
    return m_width;
}

uint32_t CmImage::getHeight() const
{
    return m_height;
}

uint32_t CmImage::getImageDataSize() const
{
    return m_imageData.size();
}

float* CmImage::getImageData()
{
    return m_imageData.data();
}

std::vector<float>& CmImage::getImageDataVector()
{
    return m_imageData;
}
//
//uint32_t CmImage::getGlTexture()
//{
//    if (m_moveToGpu)
//    {
//        m_moveToGpu = false;
//        moveToGPU_Internal();
//    }
//
//    if (CMS::usingGPU())
//    {
//        if (m_useGlobalFB && (s_framebuffer.get() != nullptr))
//        {
//            return s_framebuffer->getColorBuffer();
//        }
//        else if (!m_useGlobalFB && (m_localFramebuffer.get() != nullptr))
//        {
//            return m_localFramebuffer->getColorBuffer();
//        }
//    }
//    else if (m_texture.get())
//    {
//        return m_texture->getTexture();
//    }
//
//    return 0;
//}

void CmImage::lock()
{
    m_mutex.lock();
}

void CmImage::unlock()
{
    m_mutex.unlock();
}

void CmImage::moveToGPU()
{
    m_moveToGpu = true;
}

void CmImage::resize(uint32_t newWidth, uint32_t newHeight, bool shouldLock)
{
    if (shouldLock) lock();

    if ((m_imageData.size() > 0) && (m_width == newWidth) && (m_height == newHeight))
    {
        if (shouldLock) unlock();
        return;
    }

    if (newWidth < 1 || newHeight < 1)
    {
        if (shouldLock) unlock();
        return;
    }

    clearVector(m_imageData);

    m_width = newWidth;
    m_height = newHeight;

    m_imageData.resize(m_width * m_height * 4);

    if (shouldLock) unlock();
}

void CmImage::reset(bool shouldLock)
{
    if (shouldLock) lock();

    m_sourceName = "";

    clearVector(m_imageData);
    resize(1, 1, false);

    moveToGPU();

    if (shouldLock) unlock();
}

void CmImage::fill(std::array<float, 4> color, bool shouldLock)
{
    if (shouldLock) lock();

    for (size_t i = 0; i < m_imageData.size(); i += 4)
    {
        m_imageData[i + 0] = color[0];
        m_imageData[i + 1] = color[1];
        m_imageData[i + 2] = color[2];
        m_imageData[i + 3] = color[3];
    }

    if (shouldLock) unlock();
}

void CmImage::fill(std::vector<float> buffer, bool shouldLock)
{
    if (shouldLock) lock();

    std::copy(buffer.data(), buffer.data() + std::min(m_imageData.size(), buffer.size()), m_imageData.data());

    if (shouldLock) unlock();
}

void CmImage::fill(float* buffer, bool shouldLock)
{
    if (shouldLock) lock();

    std::copy(buffer, buffer + m_imageData.size(), m_imageData.data());

    if (shouldLock) unlock();
}

void CmImage::renderUV()
{
   // std::scoped_lock lock(m_mutex);

    int redIndex = 0;
    float u, v;
    for (size_t y = 0; y < m_height; y++)
    {
        for (size_t x = 0; x < m_width; x++)
        {
            redIndex = (y * m_width + x) * 4;
            u = ((float)x + 0.5f) / (float)m_width;
            v = ((float)y + 0.5f) / (float)m_height;
            v = 1 - v;

            m_imageData[redIndex] = u;
            m_imageData[redIndex + 1] = v;
            m_imageData[redIndex + 2] = 0;
            m_imageData[redIndex + 3] = 1;
        }
    }
}

void CmImage::moveContent(CmImage& target, bool copy)
{
   /* std::scoped_lock lock1(m_mutex);
    std::scoped_lock lock2(target);*/

    // Copy the source name
    target.m_sourceName = m_sourceName;

    // Resize the target
    target.resize(m_width, m_height, false);

    // Move / Copy the buffer
    if (copy)
        target.m_imageData = m_imageData;
    else
        target.m_imageData = std::move(m_imageData);

    // Reset self if moving
    if (!copy)
        reset(false);

    // Move to GPU
    target.moveToGPU();
}

void CmImage::moveToGPU_Internal()
{
    //std::scoped_lock lock(m_mutex);

    // Apply View Transform
    static bool lastResult = false;
    try
    {
        applyViewTransform(
            m_imageData.data(),
            m_width,
            m_height,
            m_useExposure ? CMS::getExposure() : 0.0f,
            !lastResult,
            false,
            false);
        lastResult = true;
    }
    catch (const std::runtime_error& e)
    {
        printError(__FUNCTION__, "", e.what());
        lastResult = false;
    }
}

void CmImage::applyViewTransform(
    float* buffer,
    uint32_t width,
    uint32_t height,
    float exposure,
   //std::shared_ptr<GlTexture>& texture,
   // std::shared_ptr<GlFramebuffer>& framebuffer,
    bool recreate,
    bool uploadToGPU,
    bool readback,
    std::vector<float>* outBuffer)
{
    try
    {
        CMS::ensureOK();
    }
    catch (const std::runtime_error& e)
    {
        throw std::runtime_error(makeError(__FUNCTION__, "", e.what()).c_str());
    }

    uint32_t size = width * height * 4;
    float expMul = getExposureMul(exposure);
    bool gpuMode = CMS::usingGPU() && uploadToGPU;

    // Temporary buffer to apply transforms on (CPU)
    std::vector<float> transBuffer;

    // Color Transform (CPU)
    if (!gpuMode)
    {
        transBuffer.resize(size);
        std::copy(buffer, buffer + size, transBuffer.data());

        if (exposure != 0.0f)
            for (uint32_t i = 0; i < size; i++)
                if (i % 4 != 3) transBuffer[i] *= expMul;

        try
        {
            OCIO::PackedImageDesc img(
                transBuffer.data(),
                width,
                height,
                OCIO::ChannelOrdering::CHANNEL_ORDERING_RGBA,
                OCIO::BitDepth::BIT_DEPTH_F32,
                4,                 // 4 bytes to go to the next color channel
                4 * 4,             // 4 color channels * 4 bytes per channel (till the next pixel)
                width * 4 * 4);  // width * 4 channels * 4 bytes (till the next row)

            CMS::getCpuProcessor()->apply(img);
        }
        catch (std::runtime_error& e)
        {
            throw std::runtime_error(makeError(__FUNCTION__, "Color Transform (CPU)", e.what()).c_str());
        }
    }



    // Color Transform (GPU)
  

    // Read back the result

    if (!readback || (outBuffer == nullptr)) return;

    *outBuffer = transBuffer;
      return;

}

void CmImage::cleanUp()
{
}
