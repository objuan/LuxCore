/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#if !defined(LUXCORE_DISABLE_REALBLOOM)

#include <math.h>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "luxrays/kernels/kernels.h"
#include "slg/kernels/kernels.h"
#include "slg/film/film.h"

#include "slg/film/imagepipeline/plugins/real_bloom.h"
#include "slg/film/framebuffer.h"

#include "./RealBloom/Convolution.h"
#include "./RealBloom/ColorManagement/CmImageIO.h"
#include "./RealBloom/ColorManagement/CMF.h"
#include "./RealBloom/ColorManagement/CMS.h"
#include "./RealBloom/ColorManagement/CmXYZ.h"

using namespace std;
using namespace luxrays;
using namespace slg;

using namespace RealBloom;

typedef std::unordered_map<std::string, std::string> StringMap;

struct OutputColorManagement
{
    bool applyViewTransform = false;
    std::string colorSpace = "";
    std::string display = "";
    std::string view = "";
    std::string look = "";
    float exposure = 0.0f;

    OutputColorManagement(StringMap& args, const std::string& filename) {
        colorSpace = CmImageIO::getOutputSpace();
        display = CMS::getActiveDisplay();
        view = CMS::getActiveView();
        look = CMS::getActiveLook();
    }
    void apply() {
        // Set CmImageIO parameters
        CmImageIO::setOutputSpace(colorSpace);
        CmImageIO::setApplyViewTransform(applyViewTransform);

        // Set view transform parameters
        CMS::setActiveDisplay(display);
        CMS::setActiveView(view);
        CMS::setActiveLook(look);
        CMS::setExposure(exposure);

        // Update the view transform processors if needed
        if (applyViewTransform)
            CMS::updateProcessors();
    }
};

void setInputColorSpace(const std::string& colorSpace)
{
    CmImageIO::setInputSpace(colorSpace);
    CmImageIO::setNonLinearSpace(colorSpace);
}

//------------------------------------------------------------------------------
//Intel Open Image Denoise
//------------------------------------------------------------------------------

BOOST_CLASS_EXPORT_IMPLEMENT(slg::RealBloomPlugin)

RealBloomPlugin::RealBloomPlugin(const float blendExposure, const float blendConv, const float threshold, const float knee, const std::string& knlFilename)
    :blendExposure(blendExposure), blendConv(blendConv), threshold(threshold), knee(knee), knlFilename(knlFilename), firstTime(true), imgKernel(NULL){

}

RealBloomPlugin::RealBloomPlugin() {
    blendConv =0.2f;
    blendExposure = -0.5f;
	threshold = 0.f;
	knee = 0.f;
	firstTime = true;
}

RealBloomPlugin::~RealBloomPlugin() {
    if (imgKernel != NULL)
        delete imgKernel;
}

ImagePipelinePlugin * RealBloomPlugin::Copy() const {
	return new RealBloomPlugin(blendExposure, blendConv, threshold, knee, knlFilename);
}

//------------------------------------------------------------------------------
// CPU version
//------------------------------------------------------------------------------

#if _OPENMP >= 200805
typedef unsigned int itertype;
#else
	// Visual C++ 2013 supports only OpenMP 2.5
typedef int itertype;
#endif


//bool firstTime = true;
//CmImage imgKernel;

void RealBloomPlugin::Apply(Film &film, const u_int index) {
	const double totalStartTime = WallClockTime();

    if (firstTime)
    {
        firstTime = false;

        // Load config
        //Config::load();

        // Color Management System
        if (!CMS::init())
        {
            SLG_LOG("RealBloom error");

        }

        // Color Matching Functions
        CMF::init();

        // XYZ Utility
        CmXYZ::init();

        // Color Managed Image IO
        CmImageIO::init();

        //std::string knlFilename = "C:\\Lavoro\\luxcorerender\\realbloom-main\\data\\pentagonal.exr";
        if (!boost::filesystem::exists(knlFilename))
        {
            SLG_LOG("[RealBloom] Kernel file not found " << knlFilename);

            return;
        }

        imgKernel = new CmImage();
        CmImageIO::readImage(*imgKernel, knlFilename);

        if (false)
        {
            MyConvolutionParams params;

            std::string in = "C:\\Lavoro\\luxcorerender\\realbloom-main\\data\\int04.exr";
            std::string k = "C:\\Lavoro\\luxcorerender\\realbloom-main\\data\\pentagonal.exr";

            Conv_Files(in,
                k,
                "C:\\Lavoro\\luxcorerender\\realbloom-main\\data\\out.exr",
                params, true);
        }

    }
    if (!imgKernel) return;

  //  return;

	SLG_LOG("[RealBloom] Applying ");

    Spectrum *pixels = (Spectrum *)film.channel_IMAGEPIPELINEs[index]->GetPixels();

    const u_int width = film.GetWidth();
    const u_int height = film.GetHeight();
    const u_int pixelCount = width * height;

    vector<float> rgbaBuffer(4 * pixelCount);
    vector<float> outBuffer(4 * pixelCount);

    for (u_int i = 0; i < pixelCount; ++i)
    {
        rgbaBuffer[i * 4] = pixels[i].c[0];
        rgbaBuffer[i * 4 + 1] = pixels[i].c[1];
        rgbaBuffer[i * 4 + 2] = pixels[i].c[2];
        rgbaBuffer[i * 4 + 3] = 1;
    }

    MyConvolutionParams params;
    params.blendExposure = 1;

    Conv(width, height, &rgbaBuffer,
        (int)imgKernel->getWidth(), (int)imgKernel->getHeight(), &imgKernel->getImageDataVector(),
        outBuffer.data(), "out.exr" , params, false);

    for (u_int i = 0; i < pixelCount; ++i) {
    	const u_int i4 = i * 4;
        pixels[i].c[0] = outBuffer[i4];
    	pixels[i].c[1] = outBuffer[i4+1];
    	pixels[i].c[2] = outBuffer[i4+2];
    }

	SLG_LOG("RealBloom single execution took a total of " << (boost::format("%.3f") % (WallClockTime() - totalStartTime)) << "secs");
}

void RealBloomPlugin::Conv_Files(const std::string& inpFilename, const std::string& knlFilename,
    const std::string& outFilename,
    MyConvolutionParams myParams, bool verbose) {

    CmImage imgInput;
    CmImageIO::readImage(imgInput, inpFilename);

    CmImage imgKernel;
    CmImageIO::readImage(imgKernel, knlFilename);

    float* outBuffer = new float[imgInput.getWidth() * imgInput.getHeight() * 4];

    Conv((int)imgInput.getWidth(), (int)imgInput.getHeight(), &imgInput.getImageDataVector(),
        (int)imgKernel.getWidth(), (int)imgKernel.getHeight(), &imgKernel.getImageDataVector(),
        outBuffer, outFilename,
        myParams, verbose);

    // Write the output image
    CmImage imgOut("out", "", imgInput.getWidth(), imgInput.getHeight());
    imgOut.fill(outBuffer, true);

    CmImageIO::writeImage(imgOut, outFilename);

    delete[] outBuffer;
}

void RealBloomPlugin::Conv(int in_w, int in_h, std::vector<float>* input,
    int k_w, int k_h, std::vector<float>* kernel,
    float* outBuffer, const std::string& outFilename,
    MyConvolutionParams myParams, bool verbose)
{
    //CliStackTimer totalTimer("", true);

    // Arguments

   // std::string inpFilename = args["--input"];
    std::string inpColorSpace = CMS::resolveColorSpace("Linear");// args["--input-space"]);

    // s  std::string knlFilename = args["--kernel"];
    std::string knlColorSpace = CMS::resolveColorSpace("w");// args["--kernel-space"]);

    // s  std::string outFilename = args["--output"];
    StringMap args;
    OutputColorManagement outputCM(args, outFilename);
    outputCM.applyViewTransform = true;
    outputCM.colorSpace = "Linear";
    outputCM.display = "sRGB";
    outputCM.view = "Display's Native";
    outputCM.look = "None";// AgX";
    outputCM.exposure = 0.0f;

    // Images

    CmImage imgInput("inout", "", in_w, in_h);
    CmImage imgKernel("kernel", "", k_w, k_h);
    CmImage imgConvPreview;
    CmImage imgConvResult("result", "", in_w, in_h);

    // Convolution


    RealBloom::Convolution conv;
    conv.setImgInput(&imgInput);
    conv.setImgKernel(&imgKernel);
    conv.setImgConvPreview(&imgConvPreview);
    conv.setImgConvResult(&imgConvResult);

    // Read image transform arguments
   /* readImageTransformArguments(args, "input", conv.getParams()->inputTransformParams);
    readImageTransformArguments(args, "kernel", conv.getParams()->kernelTransformParams);*/

    {
    //    CliStackTimer timer("Read the input image");
        setInputColorSpace(inpColorSpace);
        //   CmImageIO::readImage(*conv.getImgInputSrc(), inpFilename);
        conv.getImgInputSrc()->resize(in_w, in_h, true);
        conv.getImgInputSrc()->fill(*input, true);
       // timer.done(verbose);
    }

    // Read the kernel image
    {
     //   CliStackTimer timer("Read the kernel image");
        setInputColorSpace(knlColorSpace);
        conv.getImgKernelSrc()->resize(k_w, k_h, true);
        conv.getImgKernelSrc()->fill(*kernel, true);
        // CmImageIO::readImage(*conv.getImgKernelSrc(), knlFilename);
      //  timer.done(verbose);
    }

    // Parameters

    RealBloom::ConvolutionParams* params = conv.getParams();
    params->methodInfo.method = RealBloom::ConvolutionMethod::FFT_CPU;
    params->methodInfo.FFT_CPU_deconvolve = myParams.deconvolve;
    params->useKernelTransformOrigin = false,//myParams.useKernelTransformOrigin;
    params->threshold = threshold;
    params->knee = knee;
    params->autoExposure = myParams.autoExposure;
    // params->inputTransformParams = 
   //  readImageTransformArguments(args, "input", params->inputTransformParams);

    params->blendAdditive = myParams.blendAdditive;
    params->blendInput = myParams.blendInput;
    params->blendConv = blendConv;
    params->blendMix = myParams.blendMix;
    params->blendExposure = blendExposure;


    /*  std::array<float, 3> kernelColor{ 1, 1, 1 };
      float kernelRotation = 0;
      std::array<float, 2> kernelScale{ 1, 1 };
      std::array<float, 2> kernelCrop{ 1, 1 };
      std::array<float, 2> kernelCenter{ 0.5f, 0.5f };*/

    bool interrupt = false;
  
    conv.convolve();

    conv.blend();
   
    outputCM.apply();

 
    int s = in_w * in_h * 4;
       
    std::copy(imgConvResult.getImageData(), imgConvResult.getImageData() + s, outBuffer);

    //s = in_w* in_h*4;
    //for (int i = 0; i < s; i++)
    //{
    //    outBuffer[i * 4] = 1;// src[i];
    //    outBuffer[i * 4 + 1] = 1;// src[i];
    //    outBuffer[i * 4 + 2] = 0;// src[i];
    //    outBuffer[i * 4 + 3] = 1;// src[i];
    //}

}

#endif
