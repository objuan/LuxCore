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

#ifndef _SLG_INTEL_REALBLOOM_H
#define	_SLG_INTEL_REALBLOOM_H

#if !defined(LUXCORE_DISABLE_REALBLOOM)

#include <vector>
#include <string>

#include <boost/serialization/export.hpp>


#include "luxrays/luxrays.h"
#include "luxrays/core/color/color.h"
#include "slg/film/film.h"
#include "slg/film/imagepipeline/imagepipeline.h"

class CmImage;

namespace slg {

struct MyConvolutionParams
{
        /*   ImageTransformParams inputTransformParams;
        ImageTransformParams kernelTransformParams;*/
        bool useKernelTransformOrigin = true;
        float threshold = 0.0f;
        float knee = 0.0f;
        bool  autoExposure = true;
        // fisso
        bool blendAdditive = true;
        // esporre esterni
        float blendInput = 1.0f;
        // esporre esterni
        float blendConv = 1.0f;
        // esporre esterni
        float blendExposure = -5.0f;

        float blendMix = 0.2f;
        // add
        float kernelExposure = 0;
        float kernelContrast = 0;
        /* std::array<float, 3> kernelColor{ 1, 1, 1 };
         float kernelRotation = 0;
         std::array<float, 2> kernelScale{ 1, 1 };
         std::array<float, 2> kernelCrop{ 1, 1 };
         std::array<float, 2> kernelCenter{ 0.5f, 0.5f };*/
        bool  deconvolve = false;
};


//------------------------------------------------------------------------------
// Intel Open Image Denoise
//------------------------------------------------------------------------------

class RealBloomPlugin : public ImagePipelinePlugin {
public:
	RealBloomPlugin(const float blendExposure,const float blendMix, const std::string& knlFilename);

	virtual ImagePipelinePlugin *Copy() const;

	virtual void Apply(Film &film, const u_int index);

	friend class boost::serialization::access;

private:
	// Used by serialization
	RealBloomPlugin();
    ~RealBloomPlugin();

	template<class Archive> void serialize(Archive &ar, const u_int version) {
		ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ImagePipelinePlugin);
		ar & blendExposure;
        ar& blendMix;
        ar& knlFilename;
	}

	float blendExposure;
    float blendMix;
    std::string knlFilename ;
    bool firstTime ;
    CmImage* imgKernel;

	void Conv(int in_w, int in_h, std::vector<float>* input,
		int k_w, int k_h, std::vector<float>* kernel,
		float* outBuffer,
		const std::string& outFilename,
		MyConvolutionParams myParams, bool verbose);

	void Conv_Files(const std::string& inpFilename, const std::string& knlFilename,
		const std::string& outFilename,
		MyConvolutionParams myParams, bool verbose);
};

}


BOOST_CLASS_VERSION(slg::RealBloomPlugin, 1)

BOOST_CLASS_EXPORT_KEY(slg::RealBloomPlugin)

#endif
		
#endif /* _SLG_INTEL_OIDN_H */