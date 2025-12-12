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

#include "luxrays/utils/thread.h"

#include "slg/engines/pathcpu/pathcpu.h"
#include "slg/volumes/volume.h"
#include "slg/utils/varianceclamping.h"
#include "slg/samplers/metropolis.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// PathCPU RenderThread
//------------------------------------------------------------------------------

PathCPURenderThread::PathCPURenderThread(PathCPURenderEngine *engine,
		const u_int index, IntersectionDevice *device) :
		CPUNoTileRenderThread(engine, index, device) {
}

void PathCPURenderThread::RenderFunc() {
	//SLG_LOG("[PathCPURenderEngine::" << threadIndex << "] Rendering thread started");

	//--------------------------------------------------------------------------
	// Initialization
	//--------------------------------------------------------------------------

	// This is really used only by Windows for 64+ threads support
	SetThreadGroupAffinity(threadIndex);

	PathCPURenderEngine *engine = (PathCPURenderEngine *)renderEngine;
	const PathTracer &pathTracer = engine->pathTracer;

	// (engine->seedBase + 1) seed is used for sharedRndGen
	RandomGenerator *rndGen = new RandomGenerator(engine->seedBase + 1 + threadIndex);

	// Setup the sampler(s)

	Sampler *eyeSampler = nullptr;
	Sampler *lightSampler = nullptr;

	eyeSampler = engine->renderConfig->AllocSampler(rndGen, engine->film,
			nullptr, engine->samplerSharedData, Properties());
	eyeSampler->SetThreadIndex(threadIndex);
	eyeSampler->RequestSamples(PIXEL_NORMALIZED_ONLY, pathTracer.eyeSampleSize);

	if (pathTracer.hybridBackForwardEnable) {
		// Light path sampler is always Metropolis
		Properties props;
		props <<
			Property("sampler.type")("METROPOLIS") <<
			// Disable image plane meaning for samples 0 and 1
			Property("sampler.imagesamples.enable")(false) <<
			Property("sampler.metropolis.addonlycaustics")(true);

		lightSampler = Sampler::FromProperties(props, rndGen, engine->film, engine->lightSampleSplatter,
				engine->lightSamplerSharedData);
		
		lightSampler->SetThreadIndex(threadIndex);
		lightSampler->RequestSamples(SCREEN_NORMALIZED_ONLY, pathTracer.lightSampleSize);
	}

	// Setup variance clamping
	VarianceClamping varianceClamping(pathTracer.sqrtVarianceClampMaxValue);

	// Setup PathTracer thread state
	PathTracerThreadState pathTracerThreadState(device,
			eyeSampler, lightSampler,
			engine->renderConfig->scene, engine->film,
			&varianceClamping);

	// PING PONG

	PathTracerThreadState* pathTracerThreadState_1;
	RandomGenerator* ping_pong_rndGen = nullptr;
	RandomGenerator* ping_pong_rndSel = nullptr;

	if (engine->light_pingpong)
	{
		ping_pong_rndSel = new RandomGenerator(engine->seedBase + 1 + threadIndex);
		ping_pong_rndGen = new RandomGenerator(engine->seedBase + 1 + threadIndex);

		Sampler* eyeSampler_1 = nullptr;
		Sampler* lightSampler_1 = nullptr;

		eyeSampler_1 = engine->renderConfig->AllocSampler(ping_pong_rndGen, engine->film,
			nullptr, engine->samplerSharedData, Properties());
		eyeSampler_1->SetThreadIndex(threadIndex);
		eyeSampler_1->RequestSamples(PIXEL_NORMALIZED_ONLY, pathTracer.eyeSampleSize);

		if (pathTracer.hybridBackForwardEnable) {
			// Light path sampler is always Metropolis
			Properties props;
			props <<
				Property("sampler.type")("METROPOLIS") <<
				// Disable image plane meaning for samples 0 and 1
				Property("sampler.imagesamples.enable")(false) <<
				Property("sampler.metropolis.addonlycaustics")(true);

			lightSampler_1 = Sampler::FromProperties(props, ping_pong_rndGen, engine->film, engine->lightSampleSplatter,
				engine->lightSamplerSharedData);

			lightSampler_1->SetThreadIndex(threadIndex);
			lightSampler_1->RequestSamples(SCREEN_NORMALIZED_ONLY, pathTracer.lightSampleSize);
		}

		pathTracerThreadState_1 = new PathTracerThreadState(device,
			eyeSampler_1, lightSampler_1,
			engine->renderConfig->scene, engine->film,
			&varianceClamping);
	}

	// default
	if (engine->light_pingpong)
		envgroup = engine->envgroup == 0 ? LightStrategyQuery::TYPE_QUERY_GROUP_0 : LightStrategyQuery::TYPE_QUERY_GROUP_1;
	else
		envgroup = TYPE_QUERY_ALL;

	lightGroupFilter = -1;
	//Sampler* sampler = sampler_0;
	//sampler->lightGroupFilter = lightGroupFilter;
	//vector<SampleResult>& sampleResults = sampleResults_0;
	//vector<PathVertexVM>& lightPathVertices = lightPathVertices_0;
	queryMode = LightStrategyQuery::TYPE_QUERY_ALL;
	float group0Perc = (float)engine->percgroup0 / 100;

	//--------------------------------------------------------------------------
	// Trace paths
	//--------------------------------------------------------------------------

	for (u_int steps = 0; !boost::this_thread::interruption_requested(); ++steps) {
		// Check if we are in pause mode
		if (engine->pauseMode) {
			// Check every 100ms if I have to continue the rendering
			while (!boost::this_thread::interruption_requested() && engine->pauseMode)
				boost::this_thread::sleep(boost::posix_time::millisec(100));

			if (boost::this_thread::interruption_requested())
				break;
		}


		if (engine->light_pingpong)
		{
			float p01 = ping_pong_rndSel->floatValue();
			if (p01 < group0Perc)
				queryMode = LightStrategyQuery::TYPE_QUERY_GROUP_0;
			else
				queryMode = LightStrategyQuery::TYPE_QUERY_GROUP_1;

			//queryMode = LightStrategyQuery::TYPE_QUERY_GROUP_0;

			lightGroupFilter = (queryMode == LightStrategyQuery::TYPE_QUERY_GROUP_0) ? 0 : 1;
		
			PathTracerThreadState &state = (queryMode == LightStrategyQuery::TYPE_QUERY_GROUP_0) ? pathTracerThreadState : *pathTracerThreadState_1;
			
			if (state.lightSampler)
			{
				state.lightSampler->lightGroupFilter = lightGroupFilter;
				state.lightSampler[0].lightGroupFilter = lightGroupFilter;
			}
			
			state.eyeSampler->lightGroupFilter = lightGroupFilter;
			state.eyeSampleResults[0].lightGroupFilter = lightGroupFilter;

			pathTracer.RenderSample(state);
		}
		else
			pathTracer.RenderSample(pathTracerThreadState);

#ifdef WIN32
		// Work around Windows bad scheduling
		renderThread->yield();
#endif

		// Check halt conditions
		if (engine->film->GetConvergence() == 1.f)
			break;
		
		if (engine->photonGICache) {
			try {
				const u_int spp = engine->film->GetTotalEyeSampleCount() / engine->film->GetPixelCount();
				engine->photonGICache->Update(threadIndex, spp);
			} catch (boost::thread_interrupted &ti) {
				// I have been interrupted, I must stop
				break;
			}
		}
	}

	delete eyeSampler;
	delete lightSampler;
	delete rndGen;

	threadDone = true;

	// This is done to stop threads pending on barrier wait
	// inside engine->photonGICache->Update(). This can happen when an
	// halt condition is satisfied.
	if (engine->photonGICache)
		engine->photonGICache->FinishUpdate(threadIndex);

	if (engine->light_pingpong)
	{
		delete pathTracerThreadState_1;
		delete ping_pong_rndGen;
		delete ping_pong_rndSel;
	}
	//SLG_LOG("[PathCPURenderEngine::" << threadIndex << "] Rendering thread halted");
}
