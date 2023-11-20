#include "slg/engines/bidircpubgl/bidircpubgl.h"

#include <openpgl/cpp/OpenPGL.h>

using namespace std;
using namespace luxrays;
using namespace slg;


typedef Point float3;
typedef Point float2;

/* NaN-safe math ops */

//inline float path_rng_1D(KernelGlobals kg,
//	uint rng_hash,
//	int sample,
//	int dimension)
//{
//#ifdef __DEBUG_CORRELATION__
//	return (float)drand48();
//#endif
//
//	if (kernel_data.integrator.sampling_pattern == SAMPLING_PATTERN_SOBOL_BURLEY) {
//		const uint index_mask = kernel_data.integrator.sobol_index_mask;
//		return sobol_burley_sample_1D(sample, dimension, rng_hash, index_mask);
//	}
//	else {
//		return tabulated_sobol_sample_1D(kg, sample, rng_hash, dimension);
//	}
//}

//
//inline float path_state_rng_1D(KernelGlobals kg,
//	ccl_private const RNGState* rng_state,
//	const int dimension)
//{
//	return path_rng_1D(
//		kg, rng_state->rng_hash, rng_state->sample, rng_state->rng_offset + dimension);
//}
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
		o[i] = min(amin[i], max(amax[i], a[i]));
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
static pgl_vec3f guiding_vec3f(const Vector& v)
{
	return openpgl::cpp::Vector3(v.x, v.y, v.z);
}

inline float3 zero_float3() { return float3(0, 0, 0); }
inline float3 one_float3() { return float3(1,1,1); }


/* Pointers to global data structures. */
//openpgl::cpp::SampleStorage* opgl_sample_data_storage = nullptr;
//openpgl::cpp::Field* opgl_guiding_field = nullptr;
//
///* Local data structures owned by the thread. */
//openpgl::cpp::PathSegmentStorage* opgl_path_segment_storage = nullptr;
//openpgl::cpp::SurfaceSamplingDistribution* opgl_surface_sampling_distribution = nullptr;
//openpgl::cpp::VolumeSamplingDistribution* opgl_volume_sampling_distribution = nullptr;

class kernel_data
{
public:
	static bool use_guiding_direct_light;
	static bool use_guiding_mis_weights;
	//static float guiding_mis_weight;
	static int pass_stride;
	static float surface_guiding_probability;
	static float guiding_directional_sampling_type;
	static float guiding_roughness_threshold;
};


class KernelGlobals
{
public:
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
};

class IntegratorType
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

class IntegratorState
{
public:
	//Guide guiding;
	IntegratorType path;
	IntegratorType shadow_path;

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



/* Records/Adds a new path segment with the current path vertex on a surface.
 * If the path is not terminated this call is usually followed by a call of
 * guiding_record_surface_bounce. */
 void guiding_record_surface_segment(
	KernelGlobals* kg,
	IntegratorState* state,
	 Point O,Vector outDir)
{
#if  PATH_GUIDING_LEVEL >= 1
	//if (!kernel_data.integrator.train_guiding) {
	//	return;
	//}

	const pgl_vec3f zero = guiding_vec3f(zero_float3());
	const pgl_vec3f one = guiding_vec3f(one_float3());

	state->path.path_segment = kg->opgl_path_segment_storage->NextSegment();
	openpgl::cpp::SetPosition(state->path.path_segment, guiding_point3f(O));
	openpgl::cpp::SetDirectionOut(state->path.path_segment, guiding_vec3f(outDir));
	openpgl::cpp::SetVolumeScatter(state->path.path_segment, false);
	openpgl::cpp::SetScatteredContribution(state->path.path_segment, zero);
	openpgl::cpp::SetDirectContribution(state->path.path_segment, zero);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, one);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0);
#endif
}

/* Records the surface scattering event at the current vertex position of the segment. */
 void guiding_record_surface_bounce(KernelGlobals* kg,
	IntegratorState* state,
	//const Ray* ray,
	const Spectrum weight,
	const float pdf,
	const float3 N,
	const float3 wo,
	const float2 roughness,
	const float eta)
{
#if  PATH_GUIDING_LEVEL >= 4
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const float min_roughness = safe_sqrtf(fminf(roughness.x, roughness.y));
	const bool is_delta = (min_roughness == 0.0f);
	const float3 weight_rgb = spectrum_to_rgb(weight);
	const float3 normal = clamp(N, neg(one_float3()), one_float3());

	//kernel_assert(state->path.path_segment != nullptr);

	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, guiding_vec3f(one_float3()));
	openpgl::cpp::SetVolumeScatter(state->path.path_segment, false);
	openpgl::cpp::SetNormal(state->path.path_segment, guiding_vec3f(normal));
	openpgl::cpp::SetDirectionIn(state->path.path_segment, guiding_vec3f(wo));
	openpgl::cpp::SetPDFDirectionIn(state->path.path_segment, pdf);
	openpgl::cpp::SetScatteringWeight(state->path.path_segment, guiding_vec3f(weight_rgb));
	openpgl::cpp::SetIsDelta(state->path.path_segment, is_delta);
	openpgl::cpp::SetEta(state->path.path_segment, eta);
	openpgl::cpp::SetRoughness(state->path.path_segment, min_roughness);
#endif
}

/* Records the emission at the current surface intersection (physical or virtual) */
 void guiding_record_surface_emission(KernelGlobals *kg,
	IntegratorState *state,
	const Spectrum Le,
	const float mis_weight)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const float3 Le_rgb = spectrum_to_rgb(Le);

	openpgl::cpp::SetDirectContribution(state->path.path_segment, guiding_vec3f(Le_rgb));
	openpgl::cpp::SetMiWeight(state->path.path_segment, mis_weight);
#endif
}

/* Record BSSRDF Interactions */

/* Records/Adds a new path segment where the vertex position is the point of entry
 * of the sub surface scattering boundary.
 * If the path is not terminated this call is usually followed by a call of
 * guiding_record_bssrdf_weight and guiding_record_bssrdf_bounce. */
 void guiding_record_bssrdf_segment(KernelGlobals *kg,
	IntegratorState *state,
	const float3 P,
	const float3 wi)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const pgl_vec3f zero = guiding_vec3f(zero_float3());
	const pgl_vec3f one = guiding_vec3f(one_float3());

	state->path.path_segment = kg->opgl_path_segment_storage->NextSegment();
	openpgl::cpp::SetPosition(state->path.path_segment, guiding_point3f(P));
	openpgl::cpp::SetDirectionOut(state->path.path_segment, guiding_vec3f(wi));
	openpgl::cpp::SetVolumeScatter(state->path.path_segment, true);
	openpgl::cpp::SetScatteredContribution(state->path.path_segment, zero);
	openpgl::cpp::SetDirectContribution(state->path.path_segment, zero);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, one);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0);
#endif
}

/* Records the transmission of the path at the point of entry while passing
 * the surface boundary. */
void guiding_record_bssrdf_weight(KernelGlobals *kg,
	IntegratorState *state,
	const Spectrum weight,
	const Spectrum albedo)
{
#if  PATH_GUIDING_LEVEL >= 1
	//if (!kernel_data.integrator.train_guiding) {
	//	return;
	//}

	/* Note albedo left out here, will be included in guiding_record_bssrdf_bounce. */
	const float3 weight_rgb = spectrum_to_rgb(safe_divide_color(weight, albedo));

	//kernel_assert(state->path.path_segment != nullptr);

	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, guiding_vec3f(zero_float3()));
	openpgl::cpp::SetScatteringWeight(state->path.path_segment, guiding_vec3f(weight_rgb));
	openpgl::cpp::SetIsDelta(state->path.path_segment, false);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0f);
	openpgl::cpp::SetRoughness(state->path.path_segment, 1.0f);
#endif
}

/* Records the direction at the point of entry the path takes when sampling the SSS contribution.
 * If not terminated this function is usually followed by a call of
 * guiding_record_volume_transmission to record the transmittance between the point of entry and
 * the point of exit. */
void guiding_record_bssrdf_bounce(KernelGlobals *kg,
	IntegratorState *state,
	const float pdf,
	const float3 N,
	const float3 wo,
	const Spectrum weight,
	const Spectrum albedo)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const float3 normal = clamp(N, neg(one_float3()), one_float3());
	const float3 weight_rgb = spectrum_to_rgb(weight * albedo);

	//kernel_assert(state->path.path_segment != nullptr);

	openpgl::cpp::SetVolumeScatter(state->path.path_segment, false);
	openpgl::cpp::SetNormal(state->path.path_segment, guiding_vec3f(normal));
	openpgl::cpp::SetDirectionIn(state->path.path_segment, guiding_vec3f(wo));
	openpgl::cpp::SetPDFDirectionIn(state->path.path_segment, pdf);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, guiding_vec3f(weight_rgb));
#endif
}

/* Record Volume Interactions */

/* Records/Adds a new path segment with the current path vertex being inside a volume.
 * If the path is not terminated this call is usually followed by a call of
 * guiding_record_volume_bounce. */
void guiding_record_volume_segment(KernelGlobals *kg,
	IntegratorState *state,
	const float3 P,
	const float3 I)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const pgl_vec3f zero = guiding_vec3f(zero_float3());
	const pgl_vec3f one = guiding_vec3f(one_float3());

	state->path.path_segment = kg->opgl_path_segment_storage->NextSegment();

	openpgl::cpp::SetPosition(state->path.path_segment, guiding_point3f(P));
	openpgl::cpp::SetDirectionOut(state->path.path_segment, guiding_vec3f(I));
	openpgl::cpp::SetVolumeScatter(state->path.path_segment, true);
	openpgl::cpp::SetScatteredContribution(state->path.path_segment, zero);
	openpgl::cpp::SetDirectContribution(state->path.path_segment, zero);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, one);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0);
#endif
}


/* Records the volume scattering event at the current vertex position of the segment. */
void guiding_record_volume_bounce(KernelGlobals *kg,
	IntegratorState *state,
	//ccl_private const ShaderData* sd,
	const Spectrum weight,
	const float pdf,
	const float3 wo,
	const float roughness)
{
#if  PATH_GUIDING_LEVEL >= 4
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const float3 weight_rgb = spectrum_to_rgb(weight);
	const float3 normal = float3(0.0f, 0.0f, 1.0f);

	//kernel_assert(state->path.path_segment != nullptr);

	openpgl::cpp::SetVolumeScatter(state->path.path_segment, true);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, guiding_vec3f(one_float3()));
	openpgl::cpp::SetNormal(state->path.path_segment, guiding_vec3f(normal));
	openpgl::cpp::SetDirectionIn(state->path.path_segment, guiding_vec3f(wo));
	openpgl::cpp::SetPDFDirectionIn(state->path.path_segment, pdf);
	openpgl::cpp::SetScatteringWeight(state->path.path_segment, guiding_vec3f(weight_rgb));
	openpgl::cpp::SetIsDelta(state->path.path_segment, false);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0f);
	openpgl::cpp::SetRoughness(state->path.path_segment, roughness);
#endif
}

/* Records the transmission (a.k.a. transmittance weight) between the current path segment
 * and the next one, when the path is inside or passes a volume. */
void guiding_record_volume_transmission(KernelGlobals *kg,
	IntegratorState *state,
	const float3 transmittance_weight)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/

	if (state->path.path_segment) {
		// TODO (sherholz): need to find a better way to avoid this check
		if ((transmittance_weight[0] < 0.0f || !std::isfinite(transmittance_weight[0]) ||
			std::isnan(transmittance_weight[0])) ||
			(transmittance_weight[1] < 0.0f || !std::isfinite(transmittance_weight[1]) ||
				std::isnan(transmittance_weight[1])) ||
			(transmittance_weight[2] < 0.0f || !std::isfinite(transmittance_weight[2]) ||
				std::isnan(transmittance_weight[2])))
		{
		}
		else {
			openpgl::cpp::SetTransmittanceWeight(state->path.path_segment,
				guiding_vec3f(transmittance_weight));
		}
	}
#endif
}

/* Records the emission of a volume at the vertex of the current path segment. */
void guiding_record_volume_emission(KernelGlobals *kg,
	IntegratorState *state,
	const Spectrum Le)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/

	if (state->path.path_segment) {
		const float3 Le_rgb = spectrum_to_rgb(Le);

		openpgl::cpp::SetDirectContribution(state->path.path_segment, guiding_vec3f(Le_rgb));
		openpgl::cpp::SetMiWeight(state->path.path_segment, 1.0f);
	}
#endif
}

//#  define INTEGRATOR_STATE(state, nested_struct, member) \
//    kernel_integrator_state.nested_struct.member[state]

/* Record Light Interactions */

/* Adds a pseudo path vertex/segment when intersecting a virtual light source.
 * (e.g., area, sphere, or disk light). This call is often followed
 * a call of guiding_record_surface_emission, if the intersected light source
 * emits light in the direction of the path. */
void guiding_record_light_surface_segment(
	KernelGlobals *kg, IntegratorState *state)//,  const Intersection* ccl_restrict isect)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	const pgl_vec3f zero = guiding_vec3f(zero_float3());
	const pgl_vec3f one = guiding_vec3f(one_float3());

	float3 P;
	float3 ray_D;
	state->getHitRay(P, ray_D);

	/*const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
	const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
	const float3 P = ray_P + isect->t * ray_D;*/

	state->path.path_segment = kg->opgl_path_segment_storage->NextSegment();
	openpgl::cpp::SetPosition(state->path.path_segment, guiding_point3f(P));
	openpgl::cpp::SetDirectionOut(state->path.path_segment, guiding_vec3f(neg(ray_D)));
	openpgl::cpp::SetNormal(state->path.path_segment, guiding_vec3f(neg(ray_D)));
	openpgl::cpp::SetDirectionIn(state->path.path_segment, guiding_vec3f(ray_D));
	openpgl::cpp::SetPDFDirectionIn(state->path.path_segment, 1.0f);
	openpgl::cpp::SetVolumeScatter(state->path.path_segment, false);
	openpgl::cpp::SetScatteredContribution(state->path.path_segment, zero);
	openpgl::cpp::SetDirectContribution(state->path.path_segment, zero);
	openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, one);
	openpgl::cpp::SetScatteringWeight(state->path.path_segment, one);
	openpgl::cpp::SetEta(state->path.path_segment, 1.0f);
#endif
}

/* Records/Adds a final path segment when the path leaves the scene and
 * intersects with a background light (e.g., background color,
 * distant light, or env map). The vertex for this segment is placed along
 * the current ray far out the scene. */
void guiding_record_background(KernelGlobals *kg,
	IntegratorState *state,
	const Spectrum L,
	const float mis_weight)
{
#if  PATH_GUIDING_LEVEL >= 1
	//if (!kernel_data.integrator.train_guiding) {
	//	return;
	//}

	float3 P;
	float3 ray_D;
	state->getSourceRay(P, ray_D);


	const float3 L_rgb = spectrum_to_rgb(L);
	//const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
	//const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
	//const float3 P = ray_P + (1e6f) * ray_D;
	const float3 normal = make_float3(0.0f, 0.0f, 1.0f);

	openpgl::cpp::PathSegment background_segment;
	openpgl::cpp::SetPosition(&background_segment, guiding_vec3f(P));
	openpgl::cpp::SetNormal(&background_segment, guiding_vec3f(normal));
	openpgl::cpp::SetDirectionOut(&background_segment, guiding_vec3f(neg (ray_D)));
	openpgl::cpp::SetDirectContribution(&background_segment, guiding_vec3f(L_rgb));
	openpgl::cpp::SetMiWeight(&background_segment, mis_weight);
	kg->opgl_path_segment_storage->AddSegment(background_segment);
#endif
}

/* Records direct lighting from either next event estimation or a dedicated BSDF
 * sampled shadow ray. */
void guiding_record_direct_light(KernelGlobals *kg,
  IntegratorState* state)
{
#if  PATH_GUIDING_LEVEL >= 1
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/
	if (state->shadow_path.path_segment!=NULL)
	{
	/*	const Spectrum Lo = safe_divide_color(INTEGRATOR_STATE(state, shadow_path, throughput),
			INTEGRATOR_STATE(state, shadow_path, unlit_throughput));*/

		const Spectrum Lo = safe_divide_color(state->shadow_path.throughput,
			state->shadow_path.unlit_throughput);

		const float3 Lo_rgb = spectrum_to_rgb(Lo);

		const float mis_weight = state->shadow_path.guiding_mis_weight;
		//const float mis_weight = INTEGRATOR_STATE(state, shadow_path, kernel_data::guiding_mis_weight);

		if (mis_weight == 0.0f) {
			/* Scattered contribution of a next event estimation (i.e., a direct light estimate
			 * scattered at the current path vertex towards the previous vertex). */
			openpgl::cpp::AddScatteredContribution(state->shadow_path.path_segment,
				guiding_vec3f(Lo_rgb));
		}
		else {
			/* Dedicated shadow ray for BSDF sampled ray direction.
			 * The mis weight was already folded into the throughput, so need to divide it out. */
			openpgl::cpp::SetDirectContribution(state->shadow_path.path_segment,
				guiding_vec3f(Lo_rgb / mis_weight));
			openpgl::cpp::SetMiWeight(state->shadow_path.path_segment, mis_weight);
		}
	}
#endif
}

/* Record Russian Roulette */
/* Records the probability of continuing the path at the current path segment. */
void guiding_record_continuation_probability(
	KernelGlobals *kg, IntegratorState *state, const float continuation_probability)
{
#if  PATH_GUIDING_LEVEL >= 1
	//if (!kernel_data.integrator.train_guiding) {
	//	return;
	//}

	if (state->path.path_segment) {
		openpgl::cpp::SetRussianRouletteProbability(state->path.path_segment,
			continuation_probability);
	}
#endif
}

/* Path guiding debug render passes. */
#define WITH_CYCLES_DEBUG

///* Write a set of path guiding related debug information (e.g., guiding probability at first bounce) into separate rendering passes. */
void guiding_write_debug_passes(KernelGlobals *kg,
	IntegratorState *state,
	const Ray* tay)
	//float* ccl_restrict
//	render_buffer)
{
#if  PATH_GUIDING_LEVEL >= 4
#  ifdef WITH_CYCLES_DEBUG
	/*if (!kernel_data.integrator.train_guiding) {
		return;
	}*/

	
	//if (INTEGRATOR_STATE(state, path, bounce) != 0) {
//	if(state.path.bounce != 0)
//		return;
//	}
//
////	const uint32_t render_pixel_index = INTEGRATOR_STATE(state, path, render_pixel_index);
//	const uint32_t render_pixel_index = state-> path. render_pixel_index;
//
//	const uint64_t render_buffer_offset = (uint64_t)render_pixel_index *kernel_data.pass_stride;
//
//	float* buffer = render_buffer + render_buffer_offset;
//
//	if (kernel_data.film.pass_guiding_probability != PASS_UNUSED) {
//		float guiding_prob = state->path.surface_guiding_sampling_prob;
//		film_write_pass_float(buffer + kernel_data.film.pass_guiding_probability, guiding_prob);
//	}
//
//	if (kernel_data.film.pass_guiding_avg_roughness != PASS_UNUSED) {
//		float avg_roughness = 0.0f;
//		float sum_sample_weight = 0.0f;
//		for (int i = 0; i < sd->num_closure; i++) {
//			ccl_private const ShaderClosure* sc = &sd->closure[i];
//
//			if (!CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
//				continue;
//			}
//			avg_roughness += sc->sample_weight * bsdf_get_specular_roughness_squared(sc);
//			sum_sample_weight += sc->sample_weight;
//		}
//
//		avg_roughness = avg_roughness > 0.0f ? avg_roughness / sum_sample_weight : 0.0f;
//
//		film_write_pass_float(buffer + kernel_data.film.pass_guiding_avg_roughness, avg_roughness);
//	}
#  endif
#endif
}

/* Guided BSDFs */

bool guiding_bsdf_init(KernelGlobals *kg,
	IntegratorState *state,
	const float3 P,
	const float3 N,
	float& rand)
{
	// float rand_bsdf_guiding = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);

#if  PATH_GUIDING_LEVEL >= 4
	if (kg->opgl_surface_sampling_distribution->Init(
		kg->opgl_guiding_field, guiding_point3f(P), rand)) {
		kg->opgl_surface_sampling_distribution->ApplyCosineProduct(guiding_point3f(N));
		return true;
	}
#endif

	return false;
}

/// reutrn pdf
float guiding_bsdf_sample(KernelGlobals *kg,
	IntegratorState *state,
	const float2 rand_bsdf,
	float3* wo)
{
#if  PATH_GUIDING_LEVEL >= 4
	pgl_vec3f pgl_wo;
	const pgl_point2f rand = openpgl::cpp::Point2(rand_bsdf.x, rand_bsdf.y);
	const float pdf = kg->opgl_surface_sampling_distribution->SamplePDF(rand, pgl_wo);
	*wo = make_float3(pgl_wo.x, pgl_wo.y, pgl_wo.z);
	return pdf;
#else
	return 0.0f;
#endif
}

float guiding_bsdf_pdf(KernelGlobals *kg,
	IntegratorState *state,
	const float3 wo)
{
#if  PATH_GUIDING_LEVEL >= 4
	return kg->opgl_surface_sampling_distribution->PDF(guiding_vec3f(wo));
#else
	return 0.0f;
#endif
}

float guiding_surface_incoming_radiance_pdf(KernelGlobals *kg,
	IntegratorState *state,
	const float3 wo)
{
#if  PATH_GUIDING_LEVEL >= 4
	return kg->opgl_surface_sampling_distribution->IncomingRadiancePDF(guiding_vec3f(wo));
#else
	return 0.0f;
#endif
}

/* Guided Volume Phases */

bool guiding_phase_init(KernelGlobals *kg,
	IntegratorState *state,
	const float3 P,
	const float3 D,
	const float g,
	float& rand)
{
#if  PATH_GUIDING_LEVEL >= 4
	/* we do not need to guide almost delta phase functions */
	if (fabsf(g) >= 0.99f) {
		return false;
	}

	if (kg->opgl_volume_sampling_distribution->Init(
		kg->opgl_guiding_field, guiding_point3f(P), rand)) {
		kg->opgl_volume_sampling_distribution->ApplySingleLobeHenyeyGreensteinProduct(guiding_vec3f(D),
			g);
		return true;
	}
#endif

	return false;
}

float guiding_phase_sample(KernelGlobals *kg,
	IntegratorState *state,
	const float2 rand_phase,
	float3* wo)
{
#if  PATH_GUIDING_LEVEL >= 4
	pgl_vec3f pgl_wo;
	const pgl_point2f rand = openpgl::cpp::Point2(rand_phase.x, rand_phase.y);
	const float pdf = kg->opgl_volume_sampling_distribution->SamplePDF(rand, pgl_wo);
	*wo = make_float3(pgl_wo.x, pgl_wo.y, pgl_wo.z);
	return pdf;
#else
	return 0.0f;
#endif
}

float guiding_phase_pdf(KernelGlobals *kg,
	IntegratorState *state,
	const float3 wo)
{
#if  PATH_GUIDING_LEVEL >= 4
	return kg->opgl_volume_sampling_distribution->PDF(guiding_vec3f(wo));
#else
	return 0.0f;
#endif
}

// =========================================


void guiding_push_sample_data_to_global_storage(
	KernelGlobals* kg, IntegratorState* state)
{

	/* Convert the path segment representation of the random walk into radiance samples. */
#  if PATH_GUIDING_LEVEL >= 2
	const bool use_direct_light = kernel_data::use_guiding_direct_light;
	const bool use_mis_weights = kernel_data::use_guiding_mis_weights;
	kg->opgl_path_segment_storage->PrepareSamples(use_mis_weights, use_direct_light, false);
#  endif

#  if PATH_GUIDING_LEVEL >= 3
	/* Push radiance samples from current random walk/path to the global sample storage. */
	size_t num_samples = 0;
	const openpgl::cpp::SampleData* samples = kg->opgl_path_segment_storage->GetSamples(num_samples);
	kg->opgl_sample_data_storage->AddSamples(samples, num_samples);
#  endif

	/* Clear storage for the current path, to be ready for the next path. */
	kg->opgl_path_segment_storage->Clear();
}
