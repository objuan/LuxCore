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

#include "slg/engines/bidircpubgl/bidircpubgl.h"
#include "slg/engines/bidircpubgl/bidircpubglrenderstate.h"
#include "slg/samplers/sobol.h"

using namespace luxrays;
using namespace slg;
using namespace std;

//------------------------------------------------------------------------------
// BiDirCPUBGLRenderEngine
//------------------------------------------------------------------------------

BiDirCPUBGLRenderEngine::BiDirCPUBGLRenderEngine(const RenderConfig *rcfg) :
		CPUNoTileRenderEngine(rcfg), sampleSplatter(nullptr),
		photonGICache(nullptr) {
	if (rcfg->scene->camera->GetType() == Camera::STEREO)
		throw std::runtime_error("BIDIRCPU render engine doesn't support stereo camera");

	lightPathsCount = 1;
	baseRadius = 0.f;
	radiusAlpha = 0.f;

	aovWarmupSamplerSharedData = nullptr;
}

BiDirCPUBGLRenderEngine::~BiDirCPUBGLRenderEngine() {
	delete photonGICache;
	delete aovWarmupSamplerSharedData;
}

RenderState *BiDirCPUBGLRenderEngine::GetRenderState() {
	return new BiDirCPUBGLRenderState(bootStrapSeed, photonGICache);
}

void BiDirCPUBGLRenderEngine::StartLockLess() {
	const Properties &cfg = renderConfig->cfg;

	//--------------------------------------------------------------------------
	// Check to have the right sampler settings
	//--------------------------------------------------------------------------

	CheckSamplersForNoTile(RenderEngineType2String(GetType()), cfg);

	//--------------------------------------------------------------------------
	// Rendering parameters
	//--------------------------------------------------------------------------

	maxEyePathDepth = (u_int)Max(1, cfg.Get(GetDefaultProps().Get("path.maxdepth")).Get<int>());
	maxLightPathDepth = (u_int)Max(1, cfg.Get(GetDefaultProps().Get("light.maxdepth")).Get<int>());
	
	rrDepth = (u_int)Max(1, cfg.Get(GetDefaultProps().Get("path.russianroulette.depth")).Get<int>());
	rrImportanceCap = Clamp(cfg.Get(GetDefaultProps().Get("path.russianroulette.cap")).Get<float>(), 0.f, 1.f);

	// Clamping settings
	// clamping.radiance.maxvalue is the old radiance clamping, now converted in variance clamping
	sqrtVarianceClampMaxValue = cfg.Get(Property("path.clamping.radiance.maxvalue")(0.f)).Get<float>();
	if (cfg.IsDefined("path.clamping.variance.maxvalue"))
		sqrtVarianceClampMaxValue = cfg.Get(GetDefaultProps().Get("path.clamping.variance.maxvalue")).Get<float>();
	sqrtVarianceClampMaxValue = Max(0.f, sqrtVarianceClampMaxValue);

	// Albedo AOV settings
	albedoSpecularSetting = String2AlbedoSpecularSetting(cfg.Get(GetDefaultProps().Get("path.albedospecular.type")).Get<string>());
	albedoSpecularGlossinessThreshold = Max(cfg.Get(GetDefaultProps().Get("path.albedospecular.glossinessthreshold")).Get<float>(), 0.f);

	//--------------------------------------------------------------------------
	// Restore render state if there is one
	//--------------------------------------------------------------------------

	if (startRenderState) {
		// Check if the render state is of the right type
		startRenderState->CheckEngineTag(GetObjectTag());

		BiDirCPUBGLRenderState *rs = (BiDirCPUBGLRenderState *)startRenderState;

		// Use a new seed to continue the rendering
		const u_int newSeed = rs->bootStrapSeed + 1;
		SLG_LOG("Continuing the rendering with new BIDIRCPU seed: " + ToString(newSeed));
		SetSeed(newSeed);
		
		// Transfer the ownership of PhotonGI cache pointer
		photonGICache = rs->photonGICache;
		rs->photonGICache = nullptr;

		// I have to set the scene pointer in photonGICache because it is not
		// saved by serialization
		if (photonGICache)
			photonGICache->SetScene(renderConfig->scene);

		delete startRenderState;
		startRenderState = nullptr;
	}

	//--------------------------------------------------------------------------
	// Allocate PhotonGICache if enabled
	//--------------------------------------------------------------------------

	// note: photonGICache could have been restored from the render state
	if (!photonGICache) {
		photonGICache = PhotonGICache::FromProperties(renderConfig->scene, cfg);

		// photonGICache will be nullptr if the cache is disabled
		if (photonGICache)
			photonGICache->Preprocess(renderThreads.size());
	}
	
	//--------------------------------------------------------------------------
	// Albedo and Normal AOV warm up settings
	//--------------------------------------------------------------------------

	aovWarmupSPP = Max(0u, cfg.Get(GetDefaultProps().Get("path.aovs.warmup.spp")).Get<u_int>());
	if (!film->HasChannel(Film::ALBEDO) && !film->HasChannel(Film::AVG_SHADING_NORMAL))
		aovWarmupSPP = 0;
	if (aovWarmupSPP > 0)
		aovWarmupSamplerSharedData = new SobolSamplerSharedData(seedBaseGenerator.uintValue(), film);

	//--------------------------------------------------------------------------

	delete sampleSplatter;
	sampleSplatter = new FilmSampleSplatter(pixelFilter);

	CPUNoTileRenderEngine::StartLockLess();
}

void BiDirCPUBGLRenderEngine::InitFilm() {
	film->AddChannel(Film::RADIANCE_PER_PIXEL_NORMALIZED);
	film->AddChannel(Film::RADIANCE_PER_SCREEN_NORMALIZED);
	film->SetRadianceGroupCount(renderConfig->scene->lightDefs.GetLightGroupCount());
	film->SetThreadCount(renderThreads.size());
	film->Init();

	BiDirCPUBGLRenderThread::lightSampleResultsChannels.insert(Film::RADIANCE_PER_SCREEN_NORMALIZED);

	if (film->HasChannel(Film::FilmChannelType::RADIANCE_PER_PIXEL_NORMALIZED)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::RADIANCE_PER_PIXEL_NORMALIZED);
	if (film->HasChannel(Film::FilmChannelType::ALPHA)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::ALPHA);
	if (film->HasChannel(Film::FilmChannelType::DEPTH)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::DEPTH);
	if (film->HasChannel(Film::FilmChannelType::POSITION)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::POSITION);
	if (film->HasChannel(Film::FilmChannelType::GEOMETRY_NORMAL)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::GEOMETRY_NORMAL);
	if (film->HasChannel(Film::FilmChannelType::SHADING_NORMAL)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::SHADING_NORMAL);
	if (film->HasChannel(Film::FilmChannelType::MATERIAL_ID)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::MATERIAL_ID);
	if (film->HasChannel(Film::FilmChannelType::UV)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::UV);
	if (film->HasChannel(Film::FilmChannelType::OBJECT_ID)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::OBJECT_ID);
	if (film->HasChannel(Film::FilmChannelType::SAMPLECOUNT)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::SAMPLECOUNT);
	if (film->HasChannel(Film::FilmChannelType::CONVERGENCE)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::CONVERGENCE);
	if (film->HasChannel(Film::FilmChannelType::MATERIAL_ID_COLOR)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::MATERIAL_ID_COLOR);
	if (film->HasChannel(Film::FilmChannelType::ALBEDO)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::ALBEDO);
	if (film->HasChannel(Film::FilmChannelType::AVG_SHADING_NORMAL)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::AVG_SHADING_NORMAL);
	if (film->HasChannel(Film::FilmChannelType::NOISE)) BiDirCPUBGLRenderThread::eyeSampleResultsChannels.insert(Film::NOISE);

	//BiDirCPUBGLRenderThread::eyeSampleResultsChannels =
	//	Film::FilmChannels({
	//Film::RADIANCE_PER_PIXEL_NORMALIZED, Film::ALPHA, Film::DEPTH,
	//Film::POSITION, Film::GEOMETRY_NORMAL, Film::SHADING_NORMAL, Film::MATERIAL_ID,
	//Film::UV, Film::OBJECT_ID, Film::SAMPLECOUNT, Film::CONVERGENCE,
	//Film::MATERIAL_ID_COLOR, Film::ALBEDO, Film::AVG_SHADING_NORMAL, Film::NOISE
	//		});
}

void BiDirCPUBGLRenderEngine::StopLockLess() {
	CPUNoTileRenderEngine::StopLockLess();

	delete sampleSplatter;
	sampleSplatter = nullptr;

	delete photonGICache;
	photonGICache = nullptr;
}

//------------------------------------------------------------------------------
// Static methods used by RenderEngineRegistry
//------------------------------------------------------------------------------

Properties BiDirCPUBGLRenderEngine::ToProperties(const Properties &cfg) {
	return CPUNoTileRenderEngine::ToProperties(cfg) <<
			cfg.Get(GetDefaultProps().Get("renderengine.type")) <<
			cfg.Get(GetDefaultProps().Get("path.maxdepth")) <<
			cfg.Get(GetDefaultProps().Get("light.maxdepth")) <<
			cfg.Get(GetDefaultProps().Get("path.aovs.warmup.spp")) <<
			cfg.Get(GetDefaultProps().Get("path.russianroulette.depth")) <<
			cfg.Get(GetDefaultProps().Get("path.russianroulette.cap")) <<
			cfg.Get(GetDefaultProps().Get("path.clamping.variance.maxvalue")) <<
			cfg.Get(GetDefaultProps().Get("path.albedospecular.type")) <<
			cfg.Get(GetDefaultProps().Get("path.albedospecular.glossinessthreshold")) <<
			Sampler::ToProperties(cfg) <<
			PhotonGICache::ToProperties(cfg);
}

RenderEngine *BiDirCPUBGLRenderEngine::FromProperties(const RenderConfig *rcfg) {
	return new BiDirCPUBGLRenderEngine(rcfg);
}

const Properties &BiDirCPUBGLRenderEngine::GetDefaultProps() {
	static Properties props = Properties() <<
			CPUNoTileRenderEngine::GetDefaultProps() <<
			Property("renderengine.type")(GetObjectTag()) <<
			Property("path.maxdepth")(5) <<
			Property("light.maxdepth")(5) <<
			Property("path.aovs.warmup.spp")(0) <<
			Property("path.russianroulette.depth")(3) <<
			Property("path.russianroulette.cap")(.5f) <<
			Property("path.clamping.variance.maxvalue")(0.f) <<
			Property("path.albedospecular.type")("REFLECT_TRANSMIT") <<
			Property("path.albedospecular.glossinessthreshold")(.05f) <<
			PhotonGICache::GetDefaultProps();

	return props;
}
