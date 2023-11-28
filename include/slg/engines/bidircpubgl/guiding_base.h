#pragma once
#include "slg/engines/bidircpubgl/bidircpubgl.h"

#include <openpgl/cpp/OpenPGL.h>

#include <fstream>
#include <format>


using namespace std;
using namespace luxrays;
using namespace slg;

#define PATH_GUIDING_LEVEL 4

typedef Point float3;
typedef Point float2;


extern std::ofstream* log_file;

template<typename... Args> std::string string_format(const std::string& format, Args... args)
{
	int size_s = std::snprintf(nullptr, 0, format.c_str(), args...) + 1;  // Extra space for '\0'
	if (size_s <= 0) {
		throw std::runtime_error("Error during formatting.");
	}
	auto size = static_cast<size_t>(size_s);
	std::unique_ptr<char[]> buf(new char[size]);
	std::snprintf(buf.get(), size, format.c_str(), args...);
	return std::string(buf.get(), buf.get() + size - 1);  // We don't want the '\0' inside
}

inline std::string log(const char* name, float3 p)
{
	return string_format(" %s=(%f,%f,%f)", name, p.x, p.y, p.z);
}
inline std::string log2(const char* name, float2 p)
{
	return string_format(" %s=(%f,%f)", name, p.x, p.y);
}
inline std::string log(const char* name, pgl_point2f p)
{
	return string_format(" %s=(%f,%f)", name, p.x, p.y);
}
inline std::string log(const char* name, Normal p)
{
	return string_format(" %s=(%f,%f,%f)", name, p.x, p.y, p.z);
}
inline std::string log(const char* name, Vector p)
{
	return string_format(" %s=(%f,%f,%f)", name, p.x, p.y, p.z);
}
inline std::string log(const char* name, pgl_vec3f p)
{
	return string_format(" %s=(%f,%f,%f)", name, p.x, p.y, p.z);
}
inline std::string log(const char* name, float p)
{
	return string_format(" %s=%f", name, p);
}
inline std::string log(const char* name, Spectrum p)
{
	return string_format(" %s=(%f,%f,%f)", name, p.c[0], p.c[1], p.c[2] );
}

#define kernel_data (kg->data)

/* NaN-safe math ops */
//
//inline float path_rng_1D(KernelGlobals kg,
//	u_int rng_hash,
//	int sample,
//	int dimension)
//{
//#ifdef __DEBUG_CORRELATION__
//	return (float)drand48();
//#endif
//
//	if (kernel_data.integrator.sampling_pattern == SAMPLING_PATTERN_SOBOL_BURLEY) {
//		const u_int index_mask = kernel_data::sobol_index_mask;
//		return sobol_burley_sample_1D(sample, dimension, rng_hash, index_mask);
//	}
//	else {
//		return tabulated_sobol_sample_1D(kg, sample, rng_hash, dimension);
//	}
//}
//
//


inline float safe_sqrtf(float f)
{
	return sqrtf(max(f, 0.0f));
}
inline float3 spectrum_to_rgb(Spectrum s)
{
	return s.c;
}
inline float3 clamp(const float3 &a, const float3& amin, const float3& amax) {
	//return float3(a);
	float3 o;
	for (int i = 0; i < 3; i++)
	{
		o[i] = min(amax[i], max(amin[i], a[i]));
	}
	return o;
}

inline Spectrum safe_divide_color(Spectrum a, Spectrum b)
{
	for (int i = 0; i < 3; i++)
	{
		a.c[i] = b.c[i] != 0.0f ? a.c[i] / b.c[i] : 0.0f;
	}
	/*FOREACH_SPECTRUM_CHANNEL(i) {
		GET_SPECTRUM_CHANNEL(a, i) = (GET_SPECTRUM_CHANNEL(b, i) != 0.0f) ?
			GET_SPECTRUM_CHANNEL(a, i) / GET_SPECTRUM_CHANNEL(b, i) :
			0.0f;
	}*/
	return a;
}

inline float3 neg(const float3& a)
{
	return float3(-a.x, -a.y, -a.y);
}
inline float3 make_float3(float x, float y, float z) {
	return float3(x, y, z);
}

// =========================================

static pgl_vec3f guiding_vec3f(const float3 v)
{
	return openpgl::cpp::Vector3(v.x, v.y, v.z);
}

static pgl_point3f guiding_point3f(const float3 v)
{
	return openpgl::cpp::Point3(v.x, v.y, v.z);
}
static pgl_point3f guiding_point3f(const luxrays::Normal v)
{
	return openpgl::cpp::Point3(v.x, v.y, v.z);
}
static pgl_vec3f guiding_vec3f(const Vector& v)
{
	return openpgl::cpp::Vector3(v.x, v.y, v.z);
}

inline float3 zero_float3() { return float3(0, 0, 0); }
inline float3 one_float3() { return float3(1,1,1); }
inline Normal one_Normal() { return Normal(1, 1, 1); }


/* Pointers to global data structures. */
//openpgl::cpp::SampleStorage* opgl_sample_data_storage = nullptr;
//openpgl::cpp::Field* opgl_guiding_field = nullptr;
//
///* Local data structures owned by the thread. */
//openpgl::cpp::PathSegmentStorage* opgl_path_segment_storage = nullptr;
//openpgl::cpp::SurfaceSamplingDistribution* opgl_surface_sampling_distribution = nullptr;
//openpgl::cpp::VolumeSamplingDistribution* opgl_volume_sampling_distribution = nullptr;

typedef struct
{
	bool use_guiding_direct_light;
	bool use_guiding_mis_weights;
	//static float guiding_mis_weight;
	int pass_stride;
	float surface_guiding_probability=0.5;
	float guiding_directional_sampling_type;
	float guiding_roughness_threshold;

	float scrambling_distance=0;
	int seed = 0;
	
} KernelData;

typedef struct KernelGlobalsCPU 
{
	KernelData data;
	// GLOBAL
	openpgl::cpp::Field* opgl_guiding_field = nullptr;
	// generale per conservare tutti i SampleData
	openpgl::cpp::SampleStorage* opgl_sample_data_storage = nullptr;

	// FOR EACH THREAD
	
	openpgl::cpp::SurfaceSamplingDistribution* opgl_surface_sampling_distribution = nullptr;
    openpgl::cpp::VolumeSamplingDistribution* opgl_volume_sampling_distribution = nullptr;

	// utility class to help generate multiple SampleData objects during the path/random walk generation process. For the construction of a path/walk, 
	// each new PathSegment is stored in the PathSegmentStorage
	// When the walk is finished or terminated, the - radiance - SampleData is generated using a backpropagation process.The resulting samples are then be passed to the global SampleDataStorage.
	openpgl::cpp::PathSegmentStorage* opgl_path_segment_storage = nullptr;
} KernelGlobalsCPU;

typedef const KernelGlobalsCPU*  KernelGlobals;

class IntegratorPath
{
public:
	openpgl::cpp::PathSegment* path_segment;

	Spectrum throughput;
	Spectrum unlit_throughput;
	float guiding_mis_weight;

	bool bounce;
	uint32_t render_pixel_index;

	/*	const Spectrum Lo = safe_divide_color(INTEGRATOR_STATE(state, shadow_path, throughput),
		INTEGRATOR_STATE(state, shadow_path, unlit_throughput));*/
};

struct IntegratorState
{
	//Guide guiding;
	IntegratorPath path;
	IntegratorPath shadow_path;

	bool use_surface_guiding = false;
	bool use_volume_guiding = false;
	float roughness_threshold = 0.05f;
	int training_samples = 128;

	float surface_guiding_sampling_prob;
	float bssrdf_sampling_prob;
	float sample_surface_guiding_rand;

	Sampler* sampler = nullptr;

	void getSourceRay(float3& P, float3& hitDir)
	{
		/*const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
		const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
		const float3 P = ray_P + (1e6f) * ray_D;*/
	}

	void getHitRay(float3& hitP,float3 &hitDir)
	{
		/*const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
		const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
		const float3 P = ray_P + isect->t * ray_D;*/

	}
};

// =====================================
//
//float tabulated_sobol_sample_1D(KernelGlobals kg,
//	uint32_t sample,
//	const uint32_t rng_hash,
//	const uint32_t dimension)
//{
//	uint32_t seed = rng_hash;
//
//	/* Use the same sample sequence seed for all pixels when using
//   * scrambling distance. */
//	if (kernel_data.scrambling_distance < 1.0f) {
//		seed = kernel_data.seed;
//	}
//
//	/* Fetch the sample. */
//	//const uint32_t index = tabulated_sobol_shuffled_sample_index(kg, sample, dimension, seed);
//	//float x = kernel_data_fetch(sample_pattern_lut, index * NUM_TAB_SOBOL_DIMENSIONS);
//
//	///* Do limited Cranley-Patterson rotation when using scrambling distance. */
//	//if (kernel_data.integrator.scrambling_distance < 1.0f) {
//	//	const float jitter_x = hash_wang_seeded_float(dimension, rng_hash) *
//	//		kernel_data.integrator.scrambling_distance;
//	//	x += jitter_x;
//	//	x -= floorf(x);
//	//}
//
//	//return x;
//}
//
//inline float path_rng_1D(KernelGlobals kg,
//	uint32_t rng_hash,
//	int sample,
//	int dimension)
//{
//	/*if (kernel_data.integrator.sampling_pattern == SAMPLING_PATTERN_SOBOL_BURLEY) {
//		const uint32_t index_mask = kernel_data.integrator.sobol_index_mask;
//		return sobol_burley_sample_1D(sample, dimension, rng_hash, index_mask);
//	}
//	else {*/
//		return tabulated_sobol_sample_1D(kg, sample, rng_hash, dimension);
//	//}
//}

/* RNG State loaded onto stack. */
typedef struct RNGState {
	uint32_t rng_hash;
	uint32_t rng_offset;
	int sample;
} RNGState;

//inline float path_state_rng_1D(KernelGlobals kg, const RNGState* rng_state, const int dimension)
//{
//	return path_rng_1D(
//		kg, rng_state->rng_hash, rng_state->sample, rng_state->rng_offset + dimension);
//}




//class ShaderData
//{
//public:
//	float3 P;
//	Point wi;
//};

///* Records/Adds a new path segment with the current path vertex on a surface.
// * If the path is not terminated this call is usually followed by a call of * guiding_record_surface_bounce. */
//void guiding_record_surface_segment(ShaderData* sd);
//
///* Records the surface scattering event at the current vertex position of the segment. */
//void guiding_record_surface_bounce(const ShaderData* sd, const Spectrum weight, const float pdf, const float3 N, const float3 wo, const float2 roughness, const float eta);
//
///* Records the emission at the current surface intersection (physical or virtual) */
//void guiding_record_surface_emission(const Spectrum Le, const float mis_weight);
//
//
///* Records/Adds a new path segment where the vertex position is the point of entry  of the sub surface scattering boundary.
// * If the path is not terminated this call is usually followed by a call of * guiding_record_bssrdf_weight and guiding_record_bssrdf_bounce. */
//void guiding_record_bssrdf_segment(const float3 P, const float3 wi);
//
///* Records/Adds a new path segment where the vertex position is the point of entry * of the sub surface scattering boundary.
// * If the path is not terminated this call is usually followed by a call of * guiding_record_bssrdf_weight and guiding_record_bssrdf_bounce. */
//void guiding_record_bssrdf_weight(const Spectrum weight, const Spectrum albedo);
//
///* Records the direction at the point of entry the path takes when sampling the SSS contribution.
// * If not terminated this function is usually followed by a call of
// * guiding_record_volume_transmission to record the transmittance between the point of entry and  * the point of exit. */
//void guiding_record_bssrdf_bounce(const float pdf, const float3 N, const float3 wo, const Spectrum weight, const Spectrum albedo);
//
//
//
///* Records/Adds a new path segment with the current path vertex being inside a volume.
// * If the path is not terminated this call is usually followed by a call of * guiding_record_volume_bounce. */
//void guiding_record_volume_segment(const float3 P,	const float3 I);
//
///* Records the volume scattering event at the current vertex position of the segment. */
//void guiding_record_volume_bounce(const Spectrum weight, const float pdf, const float3 wo, const float roughness);
//
///* Records the direction at the point of entry the path takes when sampling the SSS contribution.
// * If not terminated this function is usually followed by a call of * guiding_record_volume_transmission to record the transmittance between the point of entry and
// * the point of exit. */
//void guiding_record_volume_transmission(const float3 transmittance_weight);
//
///* Records the emission of a volume at the vertex of the current path segment. */
//void guiding_record_volume_emission(const Spectrum Le);
//
//
//
///* Adds a pseudo path vertex/segment when intersecting a virtual light source.
// * (e.g., area, sphere, or disk light). This call is often followed * a call of guiding_record_surface_emission, if the intersected light source
// * emits light in the direction of the path. */
//void guiding_record_light_surface_segment();
//
///* Records/Adds a final path segment when the path leaves the scene and
// * intersects with a background light (e.g., background color,
// * distant light, or env map). The vertex for this segment is placed along
// * the current ray far out the scene. */
//void guiding_record_background(const Spectrum L, const float mis_weight);
//
///* Records direct lighting from either next event estimation or a dedicated BSDF
// * sampled shadow ray. */
//void guiding_record_direct_light();
//
///* Record Russian Roulette *//* Records the probability of continuing the path at the current path segment. */
//void guiding_record_continuation_probability( const float continuation_probability);
//
//
///* Guided BSDFs */
//bool  guiding_bsdf_init(const float3 P, const float3 N, float& rand);
//float guiding_bsdf_sample(const float2 rand_bsdf, float3* wo);
//float guiding_bsdf_pdf(const float3 wo);
//
//float guiding_surface_incoming_radiance_pdf(const float3 wo);
//bool guiding_phase_init(	const float3 P,	const float3 D,	const float g,	float& rand);
// float guiding_phase_pdf(const float3 wo);


