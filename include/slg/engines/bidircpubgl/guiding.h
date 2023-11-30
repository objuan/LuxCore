#pragma once
#include "slg/engines/bidircpubgl/bidircpubgl.h"
#include "guiding_base.h"


#include <openpgl/cpp/OpenPGL.h>

using namespace std;
using namespace luxrays;
using namespace slg;

/*
* global
* 
 *guiding_field,
 void *sample_data_storage

*/
// UNO PER THREAD
typedef class PathGuidingGlobalData
{
public:
	openpgl::cpp::Field* guiding_field = nullptr;
	// generale per conservare tutti i SampleData
	openpgl::cpp::SampleStorage* sample_data_storage = nullptr;

	// params
	float surface_guiding_probability = 0.5;
	int training_samples = 128;

	PathGuidingGlobalData()
	{
		openpgl::cpp::Device odevice(PGL_DEVICE_TYPE_CPU_4);

		/* The storage container which holds the training data/samples generated during the last rendering iteration. */
		sample_data_storage = new openpgl::cpp::SampleStorage();
		//guiding_field_ = make_unique<openpgl::cpp::Field>(guiding_device, field_args);

		PGLFieldArguments fieldArgs;

		pglFieldArgumentsSetDefaults(fieldArgs, PGL_SPATIAL_STRUCTURE_KDTREE, PGL_DIRECTIONAL_DISTRIBUTION_PARALLAX_AWARE_VMM);

		//	openpgl::cpp::Field* guiding_field = new openpgl::cpp::Field(odevice, field_args);

				/* The guiding field which holds the representation of the incident radiance field for the
			* complete scene. */

		guiding_field = new openpgl::cpp::Field(&odevice, fieldArgs);


	}
};

// UNO PER THREAD
typedef class KernelGlobalsCPU
{
public:
	KernelData data;

	/* The number of already performed training iterations for the guiding field. */
	int guiding_update_count = 0;

	// FROM GLOBAL
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

	KernelGlobalsCPU(PathGuidingGlobalData* globalData)
	{
		//openpgl::cpp::Device odevice(PGL_DEVICE_TYPE_CPU_4);

		///* The storage container which holds the training data/samples generated during the last rendering iteration. */
		//openpgl::cpp::SampleStorage* guiding_sample_data_storage = new openpgl::cpp::SampleStorage();
		////guiding_field_ = make_unique<openpgl::cpp::Field>(guiding_device, field_args);

		//PGLFieldArguments fieldArgs;

		//pglFieldArgumentsSetDefaults(fieldArgs, PGL_SPATIAL_STRUCTURE_KDTREE, PGL_DIRECTIONAL_DISTRIBUTION_PARALLAX_AWARE_VMM);

		////	openpgl::cpp::Field* guiding_field = new openpgl::cpp::Field(odevice, field_args);

		//		/* The guiding field which holds the representation of the incident radiance field for the
		//	* complete scene. */

		//openpgl::cpp::Field* guiding_field = new openpgl::cpp::Field(&odevice, fieldArgs);

		/* The number of already performed training iterations for the guiding field. */

		int guiding_update_count = 0;

		std::cout << "Field created successfully... exiting\n";

		//	openpgl::cpp::SurfaceSamplingDistribution *opgl_surface_sampling_distribution = new openpgl::cpp::SurfaceSamplingDistribution(guiding_field);
			//openpgl::cpp::SurfaceSamplingDistribution * opgl_volume_sampling_distribution = new openpgl::cpp::VolumeSamplingDistribution(guiding_field);

			//FIELD
			// This information can be the incidence radiance field learned from several training iterations across the whole scene
			// The Field holds separate approximations for surfaceand volumetric radiance distributions,
			// which can be accessed separately.The representation of a scene’s radiance distribution is usually separated into a positionaland directional representation
			// using a spatial subdivision structure.
			// Each spatial leaf node(a.k.a.Region) contains a directional representation for the local incident radiance distribution.

		// condivise
		opgl_sample_data_storage = globalData->sample_data_storage;
		opgl_guiding_field = globalData->guiding_field; // conserva le info della radiance

		// per thrad

		opgl_path_segment_storage = new openpgl::cpp::PathSegmentStorage();
		opgl_surface_sampling_distribution = new openpgl::cpp::SurfaceSamplingDistribution(opgl_guiding_field);
		opgl_volume_sampling_distribution = new openpgl::cpp::VolumeSamplingDistribution(opgl_guiding_field);

		int transparent_max_bounce = 128;
		int max_bounce = 256;
		bool train = true;
		if (train) {
			opgl_path_segment_storage->Reserve(transparent_max_bounce + max_bounce + 3);
			opgl_path_segment_storage->Clear();
		}
		
	}
} KernelGlobalsCPU;

typedef const KernelGlobalsCPU* KernelGlobals;

class PathGuiding
{
public:
	KernelGlobalsCPU *kg;
	IntegratorState* state;

	bool enabled;

	RNGState rng_state;

	// sampler
	Sampler* sampler = nullptr;
	int sampleStart;
	static const int SampleSize = 2;

public:

	PathGuiding(PathGuidingGlobalData* global, u_int threadIndex, Sampler* sampler, int sampleStart)
		:sampler(sampler), sampleStart(sampleStart)
	{
		//InitializeKernel();
		kg = new KernelGlobalsCPU(global);

		enabled = true;

		state = new IntegratorState();

		rng_state.rng_hash = 0;

		//int seedBase = 0; // engine->seedBaseengine->seedBase
		//RandomGenerator* rndGen = new RandomGenerator(seedBase + 1 + threadIndex);
	}

	virtual	~PathGuiding()
	{
		delete kg;
		delete state;
	}

	float path_state_rng_1D()
	{
		return sampler->GetSample(sampleStart+1);
	}

	

	void surface_shader_prepare_guiding(float3 P,Normal N)
	{
		/* Have any BSDF to guide? */
	/*	if (!(kernel_data.integrator.use_surface_guiding && (sd->flag & SD_BSDF_HAS_EVAL))) {
			state->use_surface_guiding = false;
			return;
		}*/

		const int PRNG_SURFACE_BSDF_GUIDING = 0;

		const float surface_guiding_probability = kernel_data.surface_guiding_probability;
		const int guiding_directional_sampling_type =
			kernel_data.guiding_directional_sampling_type;
		const float guiding_roughness_threshold = kernel_data.guiding_roughness_threshold;
		float rand_bsdf_guiding = sampler->GetSample(sampleStart);// path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);

		/* Compute proportion of diffuse BSDF and BSSRDFs. */
		float diffuse_sampling_fraction = 0.0f;
		float bssrdf_sampling_fraction = 0.0f;
		float bsdf_bssrdf_sampling_sum = 0.0f;

		bool fully_opaque = true;

		//for (int i = 0; i < sd->num_closure; i++) {
		//	ShaderClosure* sc = &sd->closure[i];
		//	if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
		//		const float sweight = sc->sample_weight;
		//		kernel_assert(sweight >= 0.0f);

		//		bsdf_bssrdf_sampling_sum += sweight;
		//		if (CLOSURE_IS_BSDF_DIFFUSE(sc->type) && sc->type < CLOSURE_BSDF_TRANSLUCENT_ID) {
		//			diffuse_sampling_fraction += sweight;
		//		}
		//		if (CLOSURE_IS_BSSRDF(sc->type)) {
		//			bssrdf_sampling_fraction += sweight;
		//		}

		//		if (CLOSURE_IS_BSDF_TRANSPARENT(sc->type) || CLOSURE_IS_BSDF_TRANSMISSION(sc->type)) {
		//			fully_opaque = false;
		//		}
		//	}
		//}

		//if (bsdf_bssrdf_sampling_sum > 0.0f) {
		//	diffuse_sampling_fraction /= bsdf_bssrdf_sampling_sum;
		//	bssrdf_sampling_fraction /= bsdf_bssrdf_sampling_sum;
		//}

		///* Initial guiding */
		/* The roughness because the function returns `alpha.x * alpha.y`.
		 * In addition alpha is squared again. */
		/*float avg_roughness = surface_shader_average_sample_weight_squared_roughness(sd);
		avg_roughness = safe_sqrtf(avg_roughness);
		if (!fully_opaque || avg_roughness < guiding_roughness_threshold ||
			((guiding_directional_sampling_type == GUIDING_DIRECTIONAL_SAMPLING_TYPE_PRODUCT_MIS) &&
				(diffuse_sampling_fraction <= 0.0f)) ||
			!guiding_bsdf_init(kg, state, sd->P, sd->N, rand_bsdf_guiding))
		{
			state->guiding.use_surface_guiding = false;
			state->guiding.surface_guiding_sampling_prob = 0.0f;
			return;
		}*/
		if (!fully_opaque || !guiding_bsdf_init(P, N, rand_bsdf_guiding))
		{
			state->use_surface_guiding = false;
			state->surface_guiding_sampling_prob = 0.0f;
			return;
		}

		state->use_surface_guiding = true;
		//if (kernel_data.guiding_directional_sampling_type ==
		//	GUIDING_DIRECTIONAL_SAMPLING_TYPE_PRODUCT_MIS)
		//{
		//	state->guiding.surface_guiding_sampling_prob = surface_guiding_probability *
		//		diffuse_sampling_fraction;
		//}
		//else if (kernel_data.integrator.guiding_directional_sampling_type ==
		//	GUIDING_DIRECTIONAL_SAMPLING_TYPE_RIS)
		//{
		//	state->guiding.surface_guiding_sampling_prob = surface_guiding_probability;
		//}
		//else {  // GUIDING_DIRECTIONAL_SAMPLING_TYPE_ROUGHNESS
		//	state->guiding.surface_guiding_sampling_prob = surface_guiding_probability * avg_roughness;
		//}
		state->surface_guiding_sampling_prob = surface_guiding_probability ;

		state->bssrdf_sampling_prob = bssrdf_sampling_fraction;
		state->sample_surface_guiding_rand = rand_bsdf_guiding;

		assert(state->surface_guiding_sampling_prob > 0.0f &&
			state->surface_guiding_sampling_prob <= 1.0f);
	}

	/* Records/Adds a new path segment with the current path vertex on a surface.
	 * If the path is not terminated this call is usually followed by a call of
	 * guiding_record_surface_bounce. */
	void guiding_record_surface_segment(
		Point O, Vector outDir) // world coordinate
	{
#if  PATH_GUIDING_LEVEL >= 1
		//if (!kernel_data.integrator.train_guiding) {
		//	return;
		//}

		const pgl_vec3f zero = guiding_vec3f(zero_float3());
		const pgl_vec3f one = guiding_vec3f(one_float3());

		if (log_file != NULL)
			(*log_file) << "record_surface >> add segment " << log("pos", guiding_point3f(O))
			<< log("dir", guiding_vec3f(outDir))
			<< "\n";

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
	void guiding_record_surface_bounce(
		//const Ray* ray,
		const Spectrum weight, // diviso per pdf
		const float pdf,
		const float3 N,
		const float3 wo, // world coordinate
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

		if (log_file != NULL)
			(*log_file) << "    record_surface >> bounce " << log("normal", N) << log("wo", wo)
			<< log("weight", weight)
			<< log("pdf", pdf) << log("roughness", roughness) << "\n";

		openpgl::cpp::SetTransmittanceWeight(state->path.path_segment, guiding_vec3f(one_float3()));
		openpgl::cpp::SetVolumeScatter(state->path.path_segment, false);
		openpgl::cpp::SetNormal(state->path.path_segment, guiding_vec3f(normal));
		openpgl::cpp::SetDirectionIn(state->path.path_segment, guiding_vec3f(wo));
		openpgl::cpp::SetPDFDirectionIn(state->path.path_segment, pdf);
		openpgl::cpp::SetScatteringWeight(state->path.path_segment, guiding_vec3f(weight_rgb));
		openpgl::cpp::SetIsDelta(state->path.path_segment, is_delta);
		openpgl::cpp::SetEta(state->path.path_segment, eta);
		//openpgl::cpp::SetRoughness(state->path.path_segment, min_roughness);
#endif
	}

	/* Records the emission at the current surface intersection (physical or virtual) */
	void guiding_record_surface_emission(
		const Spectrum pathThroughput,
		const Spectrum Le,
		const float mis_weight)
	{
#if  PATH_GUIDING_LEVEL >= 1
		/*if (!kernel_data.integrator.train_guiding) {
			return;
		}*/
		const float3 Le_rgb = spectrum_to_rgb(Le);

		if (log_file != NULL)
			(*log_file) << "    record_surface >> emission " << log("Le", Le)
			<< log("mis_weight", mis_weight) << "\n";

		openpgl::cpp::SetDirectContribution(state->path.path_segment, guiding_vec3f(Le_rgb));
		openpgl::cpp::SetMiWeight(state->path.path_segment, mis_weight);
#endif
	}

	/* Record BSSRDF Interactions */

	/* Records/Adds a new path segment where the vertex position is the point of entry
	 * of the sub surface scattering boundary.
	 * If the path is not terminated this call is usually followed by a call of
	 * guiding_record_bssrdf_weight and guiding_record_bssrdf_bounce. */
	void guiding_record_bssrdf_segment(
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
	void guiding_record_bssrdf_weight(
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
	void guiding_record_bssrdf_bounce(
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
	void guiding_record_volume_segment(
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
	void guiding_record_volume_bounce(
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
	void guiding_record_volume_transmission(
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
	void guiding_record_volume_emission(
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
	void guiding_record_light_surface_segment()//,  const Intersection* ccl_restrict isect)
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
	void guiding_record_background(
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

		if (log_file != NULL)
			(*log_file) << "    record_surface >> background " << log("L", L)
			<< log("mis_weight", mis_weight) << "\n";

		openpgl::cpp::PathSegment background_segment;
		openpgl::cpp::SetPosition(&background_segment, guiding_vec3f(P));
		openpgl::cpp::SetNormal(&background_segment, guiding_vec3f(normal));
		openpgl::cpp::SetDirectionOut(&background_segment, guiding_vec3f(neg(ray_D)));
		openpgl::cpp::SetDirectContribution(&background_segment, guiding_vec3f(L_rgb));
		openpgl::cpp::SetMiWeight(&background_segment, mis_weight);
		kg->opgl_path_segment_storage->AddSegment(background_segment);
#endif
	}

	/* Records direct lighting from either next event estimation or a dedicated BSDF
	 * sampled shadow ray. */
	void guiding_record_direct_light()
	{
#if  PATH_GUIDING_LEVEL >= 1
		/*if (!kernel_data.integrator.train_guiding) {
			return;
		}*/
		if (state->shadow_path.path_segment != NULL)
		{
			/*	const Spectrum Lo = safe_divide_color(INTEGRATOR_STATE(state, shadow_path, throughput),
					INTEGRATOR_STATE(state, shadow_path, unlit_throughput));*/

			const Spectrum Lo = safe_divide_color(state->shadow_path.throughput,
				state->shadow_path.unlit_throughput);

			const float3 Lo_rgb = spectrum_to_rgb(Lo);

			const float mis_weight = state->shadow_path.guiding_mis_weight;
			//const float mis_weight = INTEGRATOR_STATE(state, shadow_path, kernel_data.guiding_mis_weight);

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
	void guiding_record_continuation_probability( const float continuation_probability)
	{
#if  PATH_GUIDING_LEVEL >= 1
		//if (!kernel_data.integrator.train_guiding) {
		//	return;
		//}
		if (log_file != NULL)
			(*log_file) << log("continuation_probability", continuation_probability)
			<< "\n";

		if (state->path.path_segment) {
			openpgl::cpp::SetRussianRouletteProbability(state->path.path_segment,
				continuation_probability);
		}
#endif
	}

	/* Path guiding debug render passes. */
#define WITH_CYCLES_DEBUG

///* Write a set of path guiding related debug information (e.g., guiding probability at first bounce) into separate rendering passes. */
	void guiding_write_debug_passes(
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

	bool guiding_bsdf_init(
		const float3 P,
		const Normal N,
		float& rand)
	{
		// float rand_bsdf_guiding = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);
#if  PATH_GUIDING_LEVEL >= 4
		if (kg->opgl_surface_sampling_distribution->Init(
			kg->opgl_guiding_field, guiding_point3f(P), rand)) {
			kg->opgl_surface_sampling_distribution->ApplyCosineProduct(guiding_point3f(N));


			if (log_file != NULL)
				(*log_file) << "bsdf_init" << log("p", P) << log("N", N) << log("rand", rand) << "\n";


			return true;
		}
#endif
		if (log_file != NULL)
			(*log_file) << "bsdf_init false" << "\n";
		return false;
	}
	bool guiding_bsdf_init(
		const float3 P,
		const float3 N,
		float& rand)
	{
		// float rand_bsdf_guiding = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);

#if  PATH_GUIDING_LEVEL >= 4
		if (kg->opgl_surface_sampling_distribution->Init(
			kg->opgl_guiding_field, guiding_point3f(P), rand)) {
			kg->opgl_surface_sampling_distribution->ApplyCosineProduct(guiding_point3f(N));

			if (log_file != NULL)
				(*log_file) << "bsdf_init" << log("p", P) << log("N", N) << log("rand", rand) << "\n";
			return true;
		}
#endif
		if (log_file != NULL)
			(*log_file) << "bsdf_init false" << "\n";
		return false;
	}

//	float3 guiding_bsdf_sample_pixel(
//		const float2 point)
//	{
//#if  PATH_GUIDING_LEVEL >= 4
//		pgl_vec3f pgl_wo;
//		const pgl_point2f p = openpgl::cpp::Point2(point.x, point.y);
//		pgl_vec3f dir = kg->opgl_surface_sampling_distribution->Sample(p);
//
//		if (log_file != NULL)
//			(*log_file) << "    bsdf_sample pixel" << log("rand", rand) << "=>" << log("pdf", pdf)
//			<< log("wo", *wo)
//			<< "\n";
//
//
//		return float3(dir.x, dir.y, dir.z);
//#else
//		return 0.0f;
//#endif
//	}

	/// reutrn pdf
	float guiding_bsdf_sample(
		const float2 rand_bsdf,
		float3* wo)
	{
#if  PATH_GUIDING_LEVEL >= 4
		pgl_vec3f pgl_wo;
		const pgl_point2f rand = openpgl::cpp::Point2(rand_bsdf.x, rand_bsdf.y);
		const float pdf = kg->opgl_surface_sampling_distribution->SamplePDF(rand, pgl_wo);
		*wo = make_float3(pgl_wo.x, pgl_wo.y, pgl_wo.z);

		if (log_file != NULL)
			(*log_file) << "    bsdf_sample" << log("rand", rand) << "=>" << log("pdf", pdf)
			<< log("wo", *wo)
			<< "\n";


		return pdf;
#else
		return 0.0f;
#endif
	}


	float guiding_bsdf_pdf(const float3 wo)
	{
#if PATH_GUIDING_LEVEL >= 4
		float pdf = kg->opgl_surface_sampling_distribution->PDF(guiding_vec3f(wo));

		if (log_file != NULL)
			(*log_file) << "    bsdf_pdf" << log("wo", wo) << "=>" << log("pdf", pdf) << "\n";
		return pdf;
#else
		if (log_file != NULL)
			(*log_file) << "    bsdf_pdf" << log("wo", wo) << "=>0\n";
		return 0.0f;
#endif
	}

	float guiding_surface_incoming_radiance_pdf(
		const float3 wo)
	{
#if  PATH_GUIDING_LEVEL >= 4
		float pdf = kg->opgl_surface_sampling_distribution->IncomingRadiancePDF(guiding_vec3f(wo));

		if (log_file != NULL)
			(*log_file) << "    radiance_pdf" << log("wo", wo) << "=>" << log("pdf", pdf) << "\n";
		return pdf;
#else
		return 0.0f;
#endif
	}

	/* Guided Volume Phases */

	bool guiding_phase_init(
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

	float guiding_phase_sample(
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

	float guiding_phase_pdf(
		const float3 wo)
	{
#if  PATH_GUIDING_LEVEL >= 4
		return kg->opgl_volume_sampling_distribution->PDF(guiding_vec3f(wo));
#else
		return 0.0f;
#endif
	}

	// =========================================


	void guiding_push_sample_data_to_global_storage()
	{
		int s = kg->opgl_path_segment_storage->GetNumSegments();

		const bool validSegments = kg->opgl_path_segment_storage->ValidateSegments();
		pgl_vec3f pgl_final_color = kg->opgl_path_segment_storage->CalculatePixelEstimate(false);

		/* Convert the path segment representation of the random walk into radiance samples. */
#  if PATH_GUIDING_LEVEL >= 2
		const bool use_direct_light = kernel_data.use_guiding_direct_light;
		const bool use_mis_weights = kernel_data.use_guiding_mis_weights;
		kg->opgl_path_segment_storage->PrepareSamples(use_mis_weights, use_direct_light, false);
#  endif

#  if PATH_GUIDING_LEVEL >= 3
		/* Push radiance samples from current random walk/path to the global sample storage. */
		size_t num_samples = 0;
		const openpgl::cpp::SampleData* samples = kg->opgl_path_segment_storage->GetSamples(num_samples);
		if (num_samples > 0)
			kg->opgl_sample_data_storage->AddSamples(samples, num_samples);
#  endif

		/* Clear storage for the current path, to be ready for the next path. */
		kg->opgl_path_segment_storage->Clear();
	}
};