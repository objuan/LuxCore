#include "slg/engines/bidircpubgl/Fastfilm.h"
#include "slg/film/sampleresult.h"
#include "slg/utils/varianceclamping.h"
#include "slg/film/denoiser/filmdenoiser.h"

using namespace std;
using namespace luxrays;
using namespace slg;


namespace slg {

	//void *FillHandler()



	FastFilm::FastFilm(const u_int width, const u_int height, const u_int* subRegion) : Film(width, height, subRegion),
		filmThreadDatas(nullptr)
	{

	}
	FastFilm::~FastFilm() {
	//	~Fast();
		if (filmThreadDatas)
			delete filmThreadDatas;
	}

	void Add()
	{
	}

	void FastFilm::SetThreadCount(const u_int threadCount) {
		Film::SetThreadCount(threadCount);

		this->threadCount = threadCount;
		if (filmThreadDatas)
			delete filmThreadDatas;
		filmThreadDatas = new FilmThreadData[threadCount];
		for (int i = 0; i < threadCount; i++) {
			//filmThreadDatas[i].list = new Sample[1000];
		}
	}

	// a inizio render
	void FastFilm::Init(){
		Film::Init();

		for (int i = 0; i < threadCount; i++)
		{
			//filmThreadDatas[i].sampleCount = 0;

		/*	if (HasChannel(RADIANCE_PER_PIXEL_NORMALIZED)) {
				filmThreadDatas[i].channel_RADIANCE_PER_PIXEL_NORMALIZEDs.resize(radianceGroupCount, nullptr);
				for (u_int i = 0; i < radianceGroupCount; ++i) {
					filmThreadDatas[i].channel_RADIANCE_PER_PIXEL_NORMALIZEDs[i] = new GenericFrameBuffer<4, 1, float>(width, height);
					filmThreadDatas[i].channel_RADIANCE_PER_PIXEL_NORMALIZEDs[i]->Clear();
				}
			}
			if (HasChannel(RADIANCE_PER_SCREEN_NORMALIZED)) {
				filmThreadDatas[i].channel_RADIANCE_PER_SCREEN_NORMALIZEDs.resize(radianceGroupCount, nullptr);
				for (u_int i = 0; i < radianceGroupCount; ++i) {
					filmThreadDatas[i].channel_RADIANCE_PER_SCREEN_NORMALIZEDs[i] = new GenericFrameBuffer<3, 0, float>(width, height);
					filmThreadDatas[i].channel_RADIANCE_PER_SCREEN_NORMALIZEDs[i]->Clear();
				}
			}*/
		}
		
	}
	void FastFilm::StartFrame(const u_int threadIndex){
		//SLG_LOG("[FastFilm::" << threadIndex << "] StartFrame ");
		
		for (int i = 0; i < threadCount; i++)
		{
			//filmThreadDatas[i].sampleCount = 0;
		}
	}

	void FastFilm::EndFrame(const u_int threadIndex) {
		//SLG_LOG("[FastFilm::" << threadIndex << "] EndFrame ");

		for (int i = 0; i < threadCount; i++)
		{
		/*	FilmThreadData& data = filmThreadDatas[threadIndex];
			int c = data.sampleCount;
			for (int s = 0; s < c; s++) {
			}*/

			//SLG_LOG("[FastFilm::" << filmThreadDatas[threadIndex].sampleCount << "] cc ");

		}
	}

	void FastFilm::AddThreadSample(const u_int threadIndex, const u_int x, const u_int y,
		const SampleResult& sampleResult, const float weight) 
	{
		AddThreadSampleResultColor(threadIndex,x, y, sampleResult, weight);
		if (hasDataChannel)
			AddThreadSampleResultData(threadIndex,x, y, sampleResult);
	}


	void FastFilm::AddThreadSampleResultData(const u_int threadIndex,const u_int x, const u_int y,
		const SampleResult & sampleResult) {
	
	}

	template<class T>
	inline void AddIfValidWeighted_4_1(T* pixels, int index, const T* v,float weight)
	{
		if ( isnan(v[0]) || isinf(v[0])
			||isnan(v[1]) || isinf(v[1])
			||isnan(v[2]) || isinf(v[2])) return;
		
		T* pixel = &pixels[index];
		*pixel = *v * weight; pixel++; v++;
		*pixel = *v * weight; pixel++; v++;
		*pixel = *v * weight; pixel++; v++;
		*pixel = weight;

	}

	template<class T>
	inline void AddIfValidWeighted_3_0(T* pixels, int index, const T* v, float weight)
	{
		if (isnan(v[0]) || isinf(v[0])
			|| isnan(v[1]) || isinf(v[1])
			|| isnan(v[2]) || isinf(v[2])) return;

		T* pixel = &pixels[index];
		*pixel = *v * weight; pixel++; v++;
		*pixel = *v * weight; pixel++; v++;
		*pixel = *v * weight;

	}

	void FastFilm::AddThreadSampleResultColor(const u_int threadIndex,const u_int x, const u_int y,
		const SampleResult& sampleResult, const float weight) {
		filmDenoiser.AddSample(x, y, sampleResult, weight);

		if (isnan(weight) || isinf(weight)) return;

		//return;
		//SLG_LOG("[FastFilm::" << x << "," << y << " w:"<< weight);

		// T *pixel = &pixels[(x + y * width) * CHANNELS];
		int idx = (x + y * width);
		//int address4_1 = idx * 5;
		int address4 = idx * 4 ;
		int address3 = idx * 3;
		FilmThreadData& data = filmThreadDatas[threadIndex];

		assert(channel_RADIANCE_PER_PIXEL_NORMALIZEDs.size() == 1);
		assert(channel_RADIANCE_PER_SCREEN_NORMALIZEDs.size() == 1);

		//std::vector<GenericFrameBuffer<4, 1, float> *> channel_RADIANCE_PER_PIXEL_NORMALIZEDs;
		if (sampleResult.HasChannel(RADIANCE_PER_PIXEL_NORMALIZED))
			AddIfValidWeighted_4_1(channel_RADIANCE_PER_PIXEL_NORMALIZEDs[0]->GetPixels(), address4, sampleResult.radiance[0].c, weight);
				//channel_RADIANCE_PER_PIXEL_NORMALIZEDs[0]->AddIfValidWeightedPixel(x,y, sampleResult.radiance[0].c, weight);
				//	channel_RADIANCE_PER_PIXEL_NORMALIZEDs[0]->AddIfValidWeightedPixelFast_WEIGHT_CHANNELS(address5, sampleResult.radiance[0].c, weight);
				
		//std::vector<GenericFrameBuffer<3, 0, float>*> channel_RADIANCE_PER_SCREEN_NORMALIZEDs;
		if (sampleResult.HasChannel(RADIANCE_PER_SCREEN_NORMALIZED))
			AddIfValidWeighted_3_0(channel_RADIANCE_PER_SCREEN_NORMALIZEDs[0]->GetPixels(), address3, sampleResult.radiance[0].c, weight);
			//channel_RADIANCE_PER_SCREEN_NORMALIZEDs[0]->AddIfValidWeightedPixel(x, y, sampleResult.radiance[0].c, weight);
			//channel_RADIANCE_PER_SCREEN_NORMALIZEDs[0]->AddIfValidWeightedPixel(x,y, sampleResult.radiance[0].c, weight);
			//channel_RADIANCE_PER_SCREEN_NORMALIZEDs[0]->AddIfValidWeightedPixelFast(address3, sampleResult.radiance[0].c, weight);
		
		// Faster than HasChannel(ALPHA)
		if (channel_ALPHA && sampleResult.HasChannel(ALPHA))
			channel_ALPHA->AtomicAddIfValidWeightedPixel(x, y, &sampleResult.alpha, weight);

		if (hasComposingChannel) {
			// Faster than HasChannel(DIRECT_DIFFUSE)
			if (channel_DIRECT_DIFFUSE && sampleResult.HasChannel(DIRECT_DIFFUSE)) {
				const Spectrum c = sampleResult.directDiffuseReflect + sampleResult.directDiffuseTransmit;
				channel_DIRECT_DIFFUSE->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
			}
			if (channel_DIRECT_DIFFUSE_REFLECT && sampleResult.HasChannel(DIRECT_DIFFUSE_REFLECT))
				channel_DIRECT_DIFFUSE_REFLECT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.directDiffuseReflect.c, weight);
			if (channel_DIRECT_DIFFUSE_TRANSMIT && sampleResult.HasChannel(DIRECT_DIFFUSE_TRANSMIT))
				channel_DIRECT_DIFFUSE_TRANSMIT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.directDiffuseTransmit.c, weight);

			// Faster than HasChannel(DIRECT_GLOSSY)
			if (channel_DIRECT_GLOSSY && sampleResult.HasChannel(DIRECT_GLOSSY)) {
				const Spectrum c = sampleResult.directGlossyReflect + sampleResult.directGlossyTransmit;
				channel_DIRECT_GLOSSY->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
			}
			if (channel_DIRECT_GLOSSY_REFLECT && sampleResult.HasChannel(DIRECT_GLOSSY_REFLECT))
				channel_DIRECT_GLOSSY_REFLECT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.directGlossyReflect.c, weight);
			if (channel_DIRECT_GLOSSY_TRANSMIT && sampleResult.HasChannel(DIRECT_GLOSSY_TRANSMIT))
				channel_DIRECT_GLOSSY_TRANSMIT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.directGlossyTransmit.c, weight);

			// Faster than HasChannel(EMISSION)
			if (channel_EMISSION && sampleResult.HasChannel(EMISSION))
				channel_EMISSION->AtomicAddIfValidWeightedPixel(x, y, sampleResult.emission.c, weight);

			// Faster than HasChannel(INDIRECT_DIFFUSE)
			if (channel_INDIRECT_DIFFUSE && sampleResult.HasChannel(INDIRECT_DIFFUSE)) {
				const Spectrum c = sampleResult.indirectDiffuseReflect + sampleResult.indirectDiffuseTransmit;
				channel_INDIRECT_DIFFUSE->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
			}
			if (channel_INDIRECT_DIFFUSE_REFLECT && sampleResult.HasChannel(INDIRECT_DIFFUSE_REFLECT))
				channel_INDIRECT_DIFFUSE_REFLECT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectDiffuseReflect.c, weight);
			if (channel_INDIRECT_DIFFUSE_TRANSMIT && sampleResult.HasChannel(INDIRECT_DIFFUSE_TRANSMIT))
				channel_INDIRECT_DIFFUSE_TRANSMIT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectDiffuseTransmit.c, weight);

			// Faster than HasChannel(INDIRECT_GLOSSY)
			if (channel_INDIRECT_GLOSSY && sampleResult.HasChannel(INDIRECT_GLOSSY)) {
				const Spectrum c = sampleResult.indirectGlossyReflect + sampleResult.indirectGlossyTransmit;
				channel_INDIRECT_GLOSSY->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
			}
			if (channel_INDIRECT_GLOSSY_REFLECT && sampleResult.HasChannel(INDIRECT_GLOSSY_REFLECT))
				channel_INDIRECT_GLOSSY_REFLECT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectGlossyReflect.c, weight);
			if (channel_INDIRECT_GLOSSY_TRANSMIT && sampleResult.HasChannel(INDIRECT_GLOSSY_TRANSMIT))
				channel_INDIRECT_GLOSSY_TRANSMIT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectGlossyTransmit.c, weight);

			// Faster than HasChannel(INDIRECT_SPECULAR)
			if (channel_INDIRECT_SPECULAR && sampleResult.HasChannel(INDIRECT_SPECULAR)) {
				const Spectrum c = sampleResult.indirectSpecularReflect + sampleResult.indirectSpecularTransmit;
				channel_INDIRECT_SPECULAR->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
			}
			if (channel_INDIRECT_SPECULAR_REFLECT && sampleResult.HasChannel(INDIRECT_SPECULAR_REFLECT))
				channel_INDIRECT_SPECULAR_REFLECT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectSpecularReflect.c, weight);
			if (channel_INDIRECT_SPECULAR_TRANSMIT && sampleResult.HasChannel(INDIRECT_SPECULAR_TRANSMIT))
				channel_INDIRECT_SPECULAR_TRANSMIT->AtomicAddIfValidWeightedPixel(x, y, sampleResult.indirectSpecularTransmit.c, weight);

			// This is MATERIAL_ID_MASK and BY_MATERIAL_ID
			if (sampleResult.HasChannel(MATERIAL_ID)) {
				// MATERIAL_ID_MASK
				for (u_int i = 0; i < maskMaterialIDs.size(); ++i) {
					float pixel[2];
					pixel[0] = (sampleResult.materialID == maskMaterialIDs[i]) ? weight : 0.f;
					pixel[1] = weight;
					channel_MATERIAL_ID_MASKs[i]->AtomicAddPixel(x, y, pixel);
				}

				// BY_MATERIAL_ID
				if ((channel_RADIANCE_PER_PIXEL_NORMALIZEDs.size() > 0) && sampleResult.HasChannel(RADIANCE_PER_PIXEL_NORMALIZED)) {
					for (u_int index = 0; index < byMaterialIDs.size(); ++index) {
						Spectrum c;

						if (sampleResult.materialID == byMaterialIDs[index]) {
							// Merge all radiance groups
							for (u_int i = 0; i < Min<u_int>(sampleResult.radiance.Size(), channel_RADIANCE_PER_PIXEL_NORMALIZEDs.size()); ++i)
								c += sampleResult.radiance[i];
						}

						channel_BY_MATERIAL_IDs[index]->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
					}
				}
			}

			// Faster than HasChannel(DIRECT_SHADOW)
			if (channel_DIRECT_SHADOW_MASK && sampleResult.HasChannel(DIRECT_SHADOW_MASK))
				channel_DIRECT_SHADOW_MASK->AtomicAddIfValidWeightedPixel(x, y, &sampleResult.directShadowMask, weight);

			// Faster than HasChannel(INDIRECT_SHADOW_MASK)
			if (channel_INDIRECT_SHADOW_MASK && sampleResult.HasChannel(INDIRECT_SHADOW_MASK))
				channel_INDIRECT_SHADOW_MASK->AtomicAddIfValidWeightedPixel(x, y, &sampleResult.indirectShadowMask, weight);

			// Faster than HasChannel(IRRADIANCE)
			if (channel_IRRADIANCE && sampleResult.HasChannel(IRRADIANCE))
				channel_IRRADIANCE->AtomicAddIfValidWeightedPixel(x, y, sampleResult.irradiance.c, weight);

			// This is OBJECT_ID_MASK and BY_OBJECT_ID
			if (sampleResult.HasChannel(OBJECT_ID)) {
				// OBJECT_ID_MASK
				for (u_int i = 0; i < maskObjectIDs.size(); ++i) {
					float pixel[2];
					pixel[0] = (sampleResult.objectID == maskObjectIDs[i]) ? weight : 0.f;
					pixel[1] = weight;
					channel_OBJECT_ID_MASKs[i]->AtomicAddPixel(x, y, pixel);
				}

				// BY_OBJECT_ID
				if ((channel_RADIANCE_PER_PIXEL_NORMALIZEDs.size() > 0) && sampleResult.HasChannel(RADIANCE_PER_PIXEL_NORMALIZED)) {
					for (u_int index = 0; index < byObjectIDs.size(); ++index) {
						Spectrum c;

						if (sampleResult.objectID == byObjectIDs[index]) {
							// Merge all radiance groups
							for (u_int i = 0; i < Min<u_int>(sampleResult.radiance.Size(), channel_RADIANCE_PER_PIXEL_NORMALIZEDs.size()); ++i)
								c += sampleResult.radiance[i];
						}

						channel_BY_OBJECT_IDs[index]->AtomicAddIfValidWeightedPixel(x, y, c.c, weight);
					}
				}
			}

			// Faster than HasChannel(MATERIAL_ID_COLOR)
			if (channel_MATERIAL_ID_COLOR && sampleResult.HasChannel(MATERIAL_ID_COLOR)) {
				const u_int matID = sampleResult.materialID;
				const Spectrum matColID(
					(matID & 0x0000ffu) * (1.f / 255.f),
					((matID & 0x00ff00u) >> 8) * (1.f / 255.f),
					((matID & 0xff0000u) >> 16) * (1.f / 255.f));

				//GenericFrameBuffer<4, 1, float>* channel_MATERIAL_ID_COLOR;
				//channel_MATERIAL_ID_COLOR->AtomicAddIfValidWeightedPixel(x, y, matColID.c, weight);
				AddIfValidWeighted_4_1(channel_MATERIAL_ID_COLOR->GetPixels(), address4, matColID.c, weight);
			}

			// Faster than HasChannel(ALBEDO)
			// GenericFrameBuffer<4, 1, float> *channel_ALBEDO;
			if (channel_ALBEDO && sampleResult.HasChannel(ALBEDO))
				AddIfValidWeighted_4_1(channel_ALBEDO->GetPixels(), address4, sampleResult.albedo.c, weight);
				//channel_ALBEDO->AtomicAddIfValidWeightedPixel(x, y, sampleResult.albedo.c, weight);

			// Faster than HasChannel(AVG_SHADING_NORMAL)
			// GenericFrameBuffer<4, 1, float> *channel_AVG_SHADING_NORMAL;
			if (channel_AVG_SHADING_NORMAL && sampleResult.HasChannel(AVG_SHADING_NORMAL))
				AddIfValidWeighted_4_1(channel_AVG_SHADING_NORMAL->GetPixels(), address4, &sampleResult.shadingNormal.x, weight);
				//channel_AVG_SHADING_NORMAL->AtomicAddIfValidWeightedPixel(x, y, &sampleResult.shadingNormal.x, weight);
		}
		
	}
}

