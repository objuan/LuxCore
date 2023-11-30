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


#include "slg/engines/bidircpubgl/guiding.h"


using namespace std;
using namespace luxrays;
using namespace slg;

std::ofstream* log_file=NULL;
//PathGuiding *pathGuiding=NULL;

//bool kernel_data::use_guiding_direct_light = true;
//bool kernel_data::use_guiding_mis_weights = true;
//float kernel_data::surface_guiding_probability = 1;
Point old;

//------------------------------------------------------------------------------
// PathCPU RenderThread
//------------------------------------------------------------------------------

PathCPURenderThread::PathCPURenderThread(PathCPURenderEngine *engine,
		const u_int index, IntersectionDevice *device) :
		CPUNoTileRenderThread(engine, index, device), pathGuiding(NULL){
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

	Scene* scene = engine->renderConfig->scene;
	Camera* camera = scene->camera;

	// Setup the sampler(s)

	Sampler *eyeSampler = nullptr;
	Sampler *lightSampler = nullptr;

	eyeSampler = engine->renderConfig->AllocSampler(rndGen, engine->film,
			nullptr, engine->samplerSharedData, Properties());
	eyeSampler->SetThreadIndex(threadIndex);
	eyeSampler->RequestSamples(PIXEL_NORMALIZED_ONLY, pathTracer.eyeSampleSize+ PathGuiding::SampleSize);

	//PathGuiding::InitializeKernel();
	// ==============
	
	//log_file = new std::ofstream("C:\\Lavoro\\luxcorerender\\WindowsCompile\\Build_CMake\\LuxCore\\bin\\Debug\\guide_log.txt", std::ios_base::out );

	//	if (pathGuiding==NULL)
	float guiding_blend_factor = 0.0f;
	if (engine->renderConfig->cfg.HaveNames("path.guiding.blend_factor"))
		guiding_blend_factor = engine->renderConfig->cfg.Get("path.guiding.blend_factor").Get<float>();
	if (pathGuiding != NULL)
	{
		guiding_blend_factor = pathGuiding->kg->data.surface_guiding_probability;
		delete pathGuiding;
	}

	pathGuiding = new PathGuiding(engine->pathGuidingGlobalData, threadIndex, eyeSampler, pathTracer.eyeSampleSize);
	//engine->pathTracer.pathGuiding = pathGuiding;
	//engine->pathGuiding = pathGuiding;
	//engine->pathGuiding->enabled = guiding_blend_factor > 0.1f;
	pathGuiding->kg->data.surface_guiding_probability = guiding_blend_factor;
	pathGuiding->enabled = guiding_blend_factor > 0.1f;

	SLG_LOG("PathCPURenderThread: guiding_blend_factor" << pathGuiding->kg->data.surface_guiding_probability);
	// ==============
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

	pathTracerThreadState.pathGuiding = pathGuiding;
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

		/// <summary>
		/// ============
		/// 
		/// 
		{
			std::lock_guard<std::mutex> lock(engine->g_mutex);
			pathGuiding->guiding_push_sample_data_to_global_storage();
		}
		/// 
		if (threadIndex==0)
	//	if (engine->pathGuiding->state->use_surface_guiding)
		{
			/// 
			const size_t num_valid_samples = pathGuiding->kg->opgl_sample_data_storage->GetSizeSurface() +
				pathGuiding->kg->opgl_sample_data_storage->GetSizeVolume();
			//if (num_valid_samples >= 1024)
			if (num_valid_samples >= 10000 )
			{
				std::lock_guard<std::mutex> lock(engine->g_mutex);
				pathGuiding->kg->opgl_guiding_field->Update(*pathGuiding->kg->opgl_sample_data_storage);
				pathGuiding->kg->opgl_sample_data_storage->Clear();
			}
		}

		/// </summary>
	/*	
		if (steps > 300000)
			break;*/
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

	// ===================================

	Point probe = scene->testProbe;


	// =============
	float srand = 0;
	float3 rand_bsdf(0, 0);

	if (false)
	{
		engine->film->channel_IRRADIANCE->Clear();
		for (int x = 0; x < engine->film->GetWidth(); x++)
		{
			for (int y = 0; y < engine->film->GetHeight(); y++)
			{
				float val[] = { 0,1,0,1 };

				RayHit rayHit;
				Ray ray;

				camera->GenerateRay(0,
					x, y, &ray,
					new PathVolumeInfo(), 0, 0);

				bool hit = device->TraceRay(&ray, &rayHit);
				if (hit)
				{
					Point p;
					Normal n;

					// Check if it a triangle with bevel edges
					bool bevelContinueToTrace;
					const ExtMesh* mesh = scene->objDefs.GetSceneObject(rayHit.meshIndex)->GetExtMesh();

					Transform local2world;
					mesh->GetLocal2World(ray.time, local2world);
					n = mesh->GetGeometryNormal(local2world, rayHit.triangleIndex);


					float3 P = ray.o + ray.d * rayHit.t;
					float3 N = float3(n.x, n.y, n.z);
					float3 d = float3(ray.d.x, ray.d.y, ray.d.z);

					float rand_bsdf_guiding = 0;// eyeSampler->GetSample(2);//TIODO

					if (pathGuiding->guiding_bsdf_init(P, N, srand))
					{

						float3 dirOut;
						float guide_pdf = pathGuiding->guiding_bsdf_sample(rand_bsdf, &dirOut);

						// 
		//				float val[] = { 1,0,0,1 };
					//	float val[] = { N.x,N.y,N.z,1 };
						//float val[] = { pdf,0,0,1 };
						float val[] = { dirOut.x,dirOut.y,dirOut.z,1 };
						engine->film->channel_IRRADIANCE->SetPixel(x, y, val);
					}
					else
					{
						float val[] = { 1,0,0,1 };
						engine->film->channel_IRRADIANCE->SetPixel(x, y, val);
					}
				}
				else
				{
					float val[] = { 0,0,1,1 };
					engine->film->channel_IRRADIANCE->SetPixel(x, y, val);
				}
				//bool hit = device ? device->TraceRay(ray, rayHit) : dataSet->GetAccelerator(ACCEL_EMBREE)->Intersect(ray, rayHit);

				//guiding_bsdf_init(&kd,&state, );


			}
		}
	}

	while (false)
	{
		RayHit rayHit;
		Ray ray;
		
		if (old.x != scene->testProbe.x || old.y!= scene->testProbe.y)
		{
			old = scene->testProbe;
			camera->GenerateRay(0,
				scene->testProbe.x, scene->testProbe.y, &ray,
				new PathVolumeInfo(), 0, 0);

			bool hit = device->TraceRay(&ray, &rayHit);
			if (hit)
			{
				Point p;
				Normal n;

				// Check if it a triangle with bevel edges
				bool bevelContinueToTrace;
				const ExtMesh* mesh = scene->objDefs.GetSceneObject(rayHit.meshIndex)->GetExtMesh();

				Transform local2world;
				mesh->GetLocal2World(ray.time, local2world);
				n = mesh->GetGeometryNormal(local2world, rayHit.triangleIndex);


				float3 P = ray.o + ray.d * rayHit.t;
				float3 N = float3(n.x, n.y, n.z);
				float3 d = float3(ray.d.x, ray.d.y, ray.d.z);

				float rand_bsdf_guiding = 0;// eyeSampler->GetSample(2);//TIODO

				if (pathGuiding->guiding_bsdf_init(P, N, srand))
				{
					float3 dirOut;
					float guide_pdf = pathGuiding->guiding_bsdf_sample(rand_bsdf, &dirOut);
					SLG_LOG("PROBE: " << scene->testProbe.x << "," << scene->testProbe.y
						<< " pdf " << guide_pdf << " dir " << dirOut);
				}
			}
		}

		boost::this_thread::sleep_for(boost::chrono::milliseconds(100));
	}

	//SLG_LOG("[PathCPURenderEngine::" << threadIndex << "] Rendering thread halted");
}
