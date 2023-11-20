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

#ifndef _SLG_METROPOLIS_FAST_SAMPLER_H
#define	_SLG_METROPOLIS_FAST_SAMPLER_H

#include <string>
#include <vector>
#include <atomic>

#include "luxrays/core/randomgen.h"
#include "slg/slg.h"
#include "slg/film/film.h"
#include "slg/samplers/metropolis.h"

namespace slg {

class MetropolisFastSamplerSharedData : public MetropolisSamplerSharedData
{
};

class MetropolisFastSampler : public MetropolisSampler {
public:
	MetropolisFastSampler(luxrays::RandomGenerator *rnd, Film *film,
			const FilmSampleSplatter *flmSplatter, const bool imgSamplesEnable,
			const u_int maxRej, const float pLarge, const float imgRange,
			const bool addOnlyCstcs, MetropolisSamplerSharedData*samplerSharedData);
	virtual ~MetropolisFastSampler();

	virtual SamplerType GetType() const { return GetObjectType(); }
	virtual std::string GetTag() const { return GetObjectTag(); }
	
	virtual void NextSample(const std::vector<SampleResult>& sampleResults);

	//--------------------------------------------------------------------------
	// Static methods used by SamplerRegistry
	//--------------------------------------------------------------------------

	static SamplerType GetObjectType() { return METROPOLIS_FAST; }
	static std::string GetObjectTag() { return "METROPOLIS_FAST"; }
	static luxrays::Properties ToProperties(const luxrays::Properties& cfg);
	static Sampler* FromProperties(const luxrays::Properties& cfg, luxrays::RandomGenerator* rndGen,
		Film* film, const FilmSampleSplatter* flmSplatter, SamplerSharedData* sharedData);
	//static slg::ocl::Sampler* FromPropertiesOCL(const luxrays::Properties& cfg);
	static void AddRequiredChannels(Film::FilmChannels& channels, const luxrays::Properties& cfg);

private:
	static const luxrays::Properties& GetDefaultProps();

	void AtomicAddSampleToFilm(const u_int threadIndex,const SampleResult& sampleResult, const float weight = 1.f);
};

}

#endif	/* _SLG_METROPOLIS_SAMPLER_H */
