#ifndef _SLG_FAST_FILM_H
#define	_SLG_FAST_FILM_H

#include "slg/film/film.h"

namespace slg {


	class FastFilm : public Film {


		typedef struct Sample41 {
			GenericFrameBuffer<4, 1, float>* buffer;
			int address;
			float value;
		};

		typedef struct FilmThreadData
		{
			/*Sample* list;
			int sampleCount;*/

		/*	std::vector<GenericFrameBuffer<4, 1, float>*> channel_RADIANCE_PER_PIXEL_NORMALIZEDs;
			std::vector<GenericFrameBuffer<3, 0, float>*> channel_RADIANCE_PER_SCREEN_NORMALIZEDs;
			GenericFrameBuffer<2, 1, float>* channel_ALPHA;
			std::vector<GenericFrameBuffer<3, 0, float>*> channel_IMAGEPIPELINEs;
			GenericFrameBuffer<1, 0, float>* channel_DEPTH;
			GenericFrameBuffer<3, 0, float>* channel_POSITION;
			GenericFrameBuffer<3, 0, float>* channel_GEOMETRY_NORMAL;
			GenericFrameBuffer<3, 0, float>* channel_SHADING_NORMAL;
			GenericFrameBuffer<4, 1, float>* channel_AVG_SHADING_NORMAL;
			GenericFrameBuffer<1, 0, u_int>* channel_MATERIAL_ID;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_DIFFUSE;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_DIFFUSE_REFLECT;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_DIFFUSE_TRANSMIT;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_GLOSSY;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_GLOSSY_REFLECT;
			GenericFrameBuffer<4, 1, float>* channel_DIRECT_GLOSSY_TRANSMIT;
			GenericFrameBuffer<4, 1, float>* channel_EMISSION;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_DIFFUSE;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_DIFFUSE_REFLECT;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_DIFFUSE_TRANSMIT;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_GLOSSY;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_GLOSSY_REFLECT;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_GLOSSY_TRANSMIT;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_SPECULAR;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_SPECULAR_REFLECT;
			GenericFrameBuffer<4, 1, float>* channel_INDIRECT_SPECULAR_TRANSMIT;
			std::vector<GenericFrameBuffer<2, 1, float>*> channel_MATERIAL_ID_MASKs;
			GenericFrameBuffer<2, 1, float>* channel_DIRECT_SHADOW_MASK;
			GenericFrameBuffer<2, 1, float>* channel_INDIRECT_SHADOW_MASK;
			GenericFrameBuffer<2, 0, float>* channel_UV;
			GenericFrameBuffer<1, 0, float>* channel_RAYCOUNT;
			std::vector<GenericFrameBuffer<4, 1, float>*> channel_BY_MATERIAL_IDs;
			GenericFrameBuffer<4, 1, float>* channel_IRRADIANCE;
			GenericFrameBuffer<1, 0, u_int>* channel_OBJECT_ID;
			std::vector<GenericFrameBuffer<2, 1, float>*> channel_OBJECT_ID_MASKs;
			std::vector<GenericFrameBuffer<4, 1, float>*> channel_BY_OBJECT_IDs;
			GenericFrameBuffer<1, 0, u_int>* channel_SAMPLECOUNT;
			GenericFrameBuffer<1, 0, float>* channel_CONVERGENCE;
			GenericFrameBuffer<4, 1, float>* channel_MATERIAL_ID_COLOR;
			GenericFrameBuffer<4, 1, float>* channel_ALBEDO;
			GenericFrameBuffer<1, 0, float>* channel_NOISE;
			GenericFrameBuffer<1, 0, float>* channel_USER_IMPORTANCE;*/
		} ;

		FilmThreadData* filmThreadDatas;
		int threadCount;

	

	public:

		FastFilm(const u_int width, const u_int height, const u_int* subRegion = NULL);
		~FastFilm();

		void SetThreadCount(const u_int threadCount);
		void Init();

		void StartFrame(const u_int threadIndex);
		void EndFrame(const u_int threadIndex);
		
		void AddThreadSample(const u_int threadIndex,const u_int x, const u_int y,
			const SampleResult& sampleResult, const float weight = 1.f);

		void AddThreadSampleResultColor(const u_int threadIndex, const u_int x, const u_int y,
			const SampleResult& sampleResult, const float weight = 1.f);

		void AddThreadSampleResultData(const u_int threadIndex, const u_int x, const u_int y,
			const SampleResult& sampleResult);

		//void WriteSamples(const u_int threadIndex);

		//void AddSampleCount(const u_int threadIndex,
		//	const double RADIANCE_PER_PIXEL_NORMALIZED_count,
		//	const double RADIANCE_PER_SCREEN_NORMALIZED_count);

		//// Atomic method versions
		//void AtomicAddSample(const u_int x, const u_int y,
		//	const SampleResult& sampleResult, const float weight = 1.f);
		//void AtomicAddSampleResultColor(const u_int x, const u_int y,
		//	const SampleResult& sampleResult, const float weight);
		//void AtomicAddSampleResultData(const u_int x, const u_int y,
		//	const SampleResult& sampleResult);

	};
}
#endif