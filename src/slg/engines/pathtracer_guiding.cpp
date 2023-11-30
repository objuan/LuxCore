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

#include <boost/function.hpp>

#include "slg/engines/pathtracer.h"
#include "slg/engines/caches/photongi/photongicache.h"
#include "slg/samplers/metropolis.h"
#include "slg/utils/varianceclamping.h"

#include "slg/engines/bidircpubgl/guiding.h"

//extern PathGuiding* pathGuiding;

using namespace std;
using namespace luxrays;
using namespace slg;

static float _ShadowTerminatorAvoidanceFactor(const Normal& Ni, const Normal& Ns,
	const Vector& lightDir) {
	const float dotNsLightDir = Dot(Ns, lightDir);
	if (dotNsLightDir <= 0.f)
		return 0.f;

	const float dotNiNs = Dot(Ni, Ns);
	if (dotNiNs <= 0.f)
		return 0.f;

	const float G = Min(1.f, Dot(Ni, lightDir) / (dotNsLightDir * dotNiNs));
	if (G <= 0.f)
		return 0.f;

	const float G2 = G * G;
	const float G3 = G2 * G;

	return -G3 + G2 + G;
}

inline float reduce_add(const Spectrum a)
{
#if defined(__KERNEL_SSE__) && defined(__KERNEL_NEON__)
	__m128 t = a.m128;
	t[3] = 0.0f;
	return vaddvq_f32(t);
#else
	return (a.c[0] + a.c[1] + a.c[2]);
#endif
}

inline float average(const Spectrum a)
{
	return reduce_add(a) * (1.0f / 3.0f);
}

typedef struct BsdfEval {
	Spectrum diffuse;
	Spectrum glossy;
	Spectrum sum;
} BsdfEval;

typedef struct
{
	BSDF* bsdf;
	BSDFEvent* event;
	float3 N;
	float sample_weight = 1;
	float weight = 1;
	bool CLOSURE_IS_BSDF_OR_BSSRDF() const { return true; }
	bool CLOSURE_IS_BSDF() const { return true; }
	bool CLOSURE_IS_BSDF_DIFFUSE() const { return *event & BSDFEventType::DIFFUSE; }
	bool CLOSURE_IS_BSDF_GLOSSY() const { return *event & BSDFEventType::GLOSSY; }
	bool CLOSURE_IS_GLASS() const { return false; }

	void init(BSDF* bsdf, BSDFEvent* event) {
		this->bsdf = bsdf;
		this->event = event;
		N = float3(bsdf->hitPoint.geometryN.x, bsdf->hitPoint.geometryN.y, bsdf->hitPoint.geometryN.z);
	}
	
}ShaderClosure;

typedef struct ShaderData
{
	Vector wi;
	int num_closure;
	ShaderClosure closure[10];
}ShaderData;

// ============================================
#ifndef M_PI_F
#  define M_PI_F (3.1415926535897932f) /* pi */
#endif

inline Spectrum zero_spectrum() {
	return Spectrum
	(0, 0, 0);
}

struct GuidingRISSample {
	float3 rand;
	float2 sampled_roughness;
	float eta{ 1.0f };
	int label;
	float3 wo;
	float bsdf_pdf{ 0.0f };
	float guide_pdf{ 0.0f };
	float ris_target{ 0.0f };
	float ris_pdf{ 0.0f };
	float ris_weight{ 0.0f };

	float incoming_radiance_pdf{ 0.0f };
	BsdfEval bsdf_eval;
	float avg_bsdf_eval{ 0.0f };
	Spectrum eval = zero_spectrum();
};

inline void bsdf_eval_init(BsdfEval* eval, Spectrum value)
{
	eval->diffuse = zero_spectrum();
	eval->glossy = zero_spectrum();
	eval->sum = value;
}
inline float2 float3_to_float2(const float3 a)
{
	return float2(a.x, a.y);
}
inline float dot(const float3 a, const float3 b)
{
#  if defined(__KERNEL_SSE41__) && defined(__KERNEL_SSE__)
	return _mm_cvtss_f32(_mm_dp_ps(a, b, 0x7F));
#  else
	return a.x * b.x + a.y * b.y + a.z * b.z;
#  endif
}

inline float reduce_min(float3 a)
{
	return min(min(a.x, a.y), a.z);
}

/* --------------------------------------------------------------------
 * BSDF Evaluation
 *
 * BSDF evaluation result, split between diffuse and glossy. This is used to
 * accumulate render passes separately. Note that reflection, transmission
 * and volume scattering are written to different render passes, but we assume
 * that only one of those can happen at a bounce, and so do not need to accumulate
 * them separately. */

inline void bsdf_eval_init( BsdfEval* eval,
	const ShaderClosure* sc,
	const float3 wo,
	Spectrum value)
{
	eval->diffuse = zero_spectrum();
	eval->glossy = zero_spectrum();

	if (sc->CLOSURE_IS_BSDF_DIFFUSE()) {
		eval->diffuse = value;
	}
	else if (sc->CLOSURE_IS_BSDF_GLOSSY()) {
		eval->glossy = value;
	}
	else if (sc->CLOSURE_IS_GLASS()) {
		/* Glass can count as glossy or transmission, depending on which side we end up on. */
		if (dot(sc->N, wo) > 0.0f) {
			eval->glossy = value;
		}
	}
	eval->sum = value;
}


inline void bsdf_eval_accum(BsdfEval* eval,
	const ShaderClosure* sc,
	const float3 wo,
	Spectrum value)
{
	if (sc->CLOSURE_IS_BSDF_DIFFUSE()) {
		eval->diffuse += value;
	}
	else if (sc->CLOSURE_IS_BSDF_GLOSSY()) {
		eval->glossy += value;
	}
	else if (sc->CLOSURE_IS_GLASS()) {
		if (dot(sc->N, wo) > 0.0f) {
			eval->glossy += value;
		}
	}

	eval->sum += value;
}



// =====
// non / pdf
Spectrum bsdf_eval(KernelGlobals kg,
	ShaderData* sd,
	const ShaderClosure* sc,
	const float3 wo, // world coordinate
	float* pdf)
{
	Spectrum eval = zero_spectrum();
	*pdf = 0.f;

	float directPDF;
	float reversePDF;
	//Vector wo = Vector(_wo.x, _wo.x, _wo.x);

	/*Vector localEyeDir = sc->bsdf->GetFrame().ToLocal(Vector(sc->bsdf->hitPoint.fixedDir.x,
		sc->bsdf->hitPoint.fixedDir.y,
		sc->bsdf->hitPoint.fixedDir.z));*/
	Vector localEyeDir = sc->bsdf->GetFrame().ToLocal(sd->wi);
	Vector localLightDir = sc->bsdf->GetFrame().ToLocal(Vector(wo.x, wo.y, wo.z));

	eval = sc->bsdf->GetMaterial()->Evaluate(sc->bsdf->hitPoint, localLightDir,
		localEyeDir , sc->event, &directPDF,&reversePDF);
	*pdf = directPDF;
	return eval;// / directPDF;
}

inline float surface_shader_bsdf_eval_pdfs(const KernelGlobals kg,
	ShaderData* sd,
	const float3 wo,
	BsdfEval* result_eval,
	float* pdfs,
	const int light_shader_flags)
{
	/* This is the veach one-sample model with balance heuristic, some pdf
	 * factors drop out when using balance heuristic weighting. */
	float sum_pdf = 0.0f;
	float sum_sample_weight = 0.0f;
	bsdf_eval_init(result_eval, zero_spectrum());
	for (int i = 0; i < sd->num_closure; i++) {
		const ShaderClosure* sc = &sd->closure[i];

		if (sc->CLOSURE_IS_BSDF_OR_BSSRDF()) {
			if (sc->CLOSURE_IS_BSDF() )
				//&& !_surface_shader_exclude(sc->type, light_shader_flags)) 
			{
				float bsdf_pdf = 0.0f;
				Spectrum eval = bsdf_eval(kg, sd, sc, wo, &bsdf_pdf);
				if(bsdf_pdf < 0)
				{
					bsdf_pdf = -bsdf_pdf;
					int f = 0;
				}
				if (!eval.Black())
				{
					assert(bsdf_pdf >= 0.0f);
					if (bsdf_pdf != 0.0f) {
						bsdf_eval_accum(result_eval, sc, wo, eval * sc->weight);
						sum_pdf += bsdf_pdf * sc->sample_weight;
						assert(bsdf_pdf * sc->sample_weight >= 0.0f);
						pdfs[i] = bsdf_pdf * sc->sample_weight;
					}
					else {
						pdfs[i] = 0.0f;
					}
				}
				else {
					pdfs[i] = 0.0f;
				}
			}
			else {
				pdfs[i] = 0.0f;
			}

			sum_sample_weight += sc->sample_weight;
		}
		else {
			pdfs[i] = 0.0f;
		}
	}
	if (sum_pdf > 0.0f) {
		for (int i = 0; i < sd->num_closure; i++) {
			pdfs[i] /= sum_pdf;
		}
	}

	return (sum_sample_weight > 0.0f) ? sum_pdf / sum_sample_weight : 0.0f;
}

inline bool calculate_ris_target(GuidingRISSample* ris_sample,
	const float guiding_sampling_prob)
{
	const float pi_factor = 2.0f;
	if (ris_sample->avg_bsdf_eval > 0.0f && ris_sample->bsdf_pdf > 1e-10f &&
		ris_sample->guide_pdf > 0.0f)
	{
		ris_sample->ris_target = (ris_sample->avg_bsdf_eval *
			((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * M_PI_F))) +
				(guiding_sampling_prob * ris_sample->incoming_radiance_pdf))));
		ris_sample->ris_pdf = (0.5f * (ris_sample->bsdf_pdf + ris_sample->guide_pdf));
		ris_sample->ris_weight = ris_sample->ris_target / ris_sample->ris_pdf;
		return true;
	}
	ris_sample->ris_target = 0.0f;
	ris_sample->ris_pdf = 0.0f;
	return false;
}




bool SampleBSDF_guiding(PathGuiding *pathGuiding,
	BSDF* bsdf, //Vector* sampledDir,
	const float u0, const float u1,
//	Spectrum* out_result,
	//float* pdfW, 
	//float* absCosSampledDir,
	BSDFEvent* event,
	// --
	BsdfEval* bsdf_eval, // non diviso per pdf
	float3* wo, // world coordinate
	float* bsdf_pdf,
	float* mis_pdf,
	float* unguided_bsdf_pdf,
	float2* sampled_roughness,
	float* eta)
{
	ShaderData sd;
	sd.num_closure = 1;
	sd.closure[0].init(bsdf,event);
	sd.wi = bsdf->hitPoint.fixedDir;
	ShaderClosure* sc = &sd.closure[0];
	

	// =====================
	// 
	//const bool use_surface_guiding = state->guiding.use_surface_guiding;
	const float guiding_sampling_prob = pathGuiding->state->surface_guiding_sampling_prob;
	const float bssrdf_sampling_prob = pathGuiding->state->bssrdf_sampling_prob;
	float rand_bsdf_guiding = pathGuiding->state->sample_surface_guiding_rand;

	Spectrum eval = zero_spectrum();
	bsdf_eval_init(bsdf_eval, eval);

	*unguided_bsdf_pdf = 0.0f;
	float guide_pdf = 0.0f;

	// selected RIS candidate
	int ris_idx = 0;

	// meta data for the two RIS candidates
	GuidingRISSample ris_samples[2];
	ris_samples[0].rand = float3(u0, u1, u0);
	ris_samples[1].rand = float3(u0, u1, u0);

	// ----------------------------------------------------
	// generate the first RIS candidate using a BSDF sample
	// ---------------------------------------------------- 

	HitPoint& hitPoint = bsdf->hitPoint;

	Vector localFixedDir = bsdf->GetFrame().ToLocal(hitPoint.fixedDir);
	Vector localSampledDir;

	// bsdf_sample
	float pdfW0;
	// result0 è diviso per /pdf
	Spectrum result0 = bsdf->GetMaterial()->Sample(hitPoint,
		localFixedDir, &localSampledDir, u0, u1, hitPoint.passThroughEvent,
		&pdfW0, event);

	/*const Vector& localLightDir = hitPoint.fromLight ? localFixedDir : localSampledDir;
	const Vector& localEyeDir = hitPoint.fromLight ? localSampledDir : localFixedDir;*/

	ris_samples[0].eval = result0 * pdfW0; // lo riporto senza pdf
	// in world coordinate
	Vector wSampledDir = bsdf->GetFrame().ToWorld(localSampledDir);
	ris_samples[0].wo = float3(wSampledDir.x, wSampledDir.y, wSampledDir.z);
	ris_samples[0].bsdf_pdf = pdfW0;
	ris_samples[0].sampled_roughness = float2(1, 1);
	ris_samples[0].eta = 1;
	
	bsdf_eval_init(
		&ris_samples[0].bsdf_eval, sc, ris_samples[0].wo, ris_samples[0].eval * sc->weight);

	if (ris_samples[0].bsdf_pdf > 0.0f) {
	/*	if (sd->num_closure > 1) {
			float sweight = sc->sample_weight;
			ris_samples[0].bsdf_pdf = _surface_shader_bsdf_eval_mis(kg,
				sd,
				ris_samples[0].wo,
				sc,
				&ris_samples[0].bsdf_eval,
				(ris_samples[0].bsdf_pdf) *
				sweight,
				sweight,
				0);
			kernel_assert(reduce_min(bsdf_eval_sum(&ris_samples[0].bsdf_eval)) >= 0.0f);
		}*/

		// media di ris_samples[0].bsdf_eval.sum (x,y,z)
		ris_samples[0].avg_bsdf_eval = average(ris_samples[0].bsdf_eval.sum);
		ris_samples[0].guide_pdf = pathGuiding->guiding_bsdf_pdf(ris_samples[0].wo);
		ris_samples[0].guide_pdf *= (1.0f - bssrdf_sampling_prob);
		ris_samples[0].incoming_radiance_pdf = pathGuiding->guiding_surface_incoming_radiance_pdf(
			ris_samples[0].wo);
		ris_samples[0].bsdf_pdf = max(0.0f, ris_samples[0].bsdf_pdf);
	}

	///
	//bsdf_eval_init(bsdf_eval, ris_samples[0].eval);
	//*bsdf_pdf = ris_samples[0].bsdf_pdf;
	////*absCosSampledDir = fabsf(CosTheta(localSampledDir));
	//*wo = ris_samples[0].wo;
	////*event = event

	//return true;

	// ------------------------------------------------------------------------------
	// generate the second RIS candidate using a sample from the guiding distribution
	// ------------------------------------------------------------------------------

	float unguided_bsdf_pdfs[2];
	bsdf_eval_init(&ris_samples[1].bsdf_eval, eval);
	ris_samples[1].guide_pdf = pathGuiding->guiding_bsdf_sample(
		 float3_to_float2(ris_samples[1].rand), &ris_samples[1].wo);

	// to local
	//auto v =  bsdf->GetFrame().ToLocal(Vector(ris_samples[1].wo.x, ris_samples[1].wo.y, ris_samples[1].wo.z));
	//ris_samples[1].wo = float3(v.x,v.y,v.z);

	ris_samples[1].guide_pdf *= (1.0f - bssrdf_sampling_prob);
	ris_samples[1].incoming_radiance_pdf = pathGuiding->guiding_surface_incoming_radiance_pdf(
		ris_samples[1].wo);
	ris_samples[1].bsdf_pdf = surface_shader_bsdf_eval_pdfs(pathGuiding->kg, &sd,
		 ris_samples[1].wo, &ris_samples[1].bsdf_eval, unguided_bsdf_pdfs, 0);
	ris_samples[1].label = ris_samples[0].label;
	ris_samples[1].avg_bsdf_eval = average(ris_samples[1].bsdf_eval.sum);
	ris_samples[1].bsdf_pdf = max(0.0f, ris_samples[1].bsdf_pdf);

	// ------------------------------------------------------------------------------
	// calculate the RIS target functions for each RIS candidate
	// ------------------------------------------------------------------------------

	int num_ris_candidates = 0;
	float sum_ris_weights = 0.0f;
	if (calculate_ris_target(&ris_samples[0], guiding_sampling_prob)) {
		sum_ris_weights += ris_samples[0].ris_weight;
		num_ris_candidates++;
	}
	assert(ris_samples[0].ris_weight >= 0.0f);
	assert(sum_ris_weights >= 0.0f);

	if (calculate_ris_target(&ris_samples[1], guiding_sampling_prob)) {
		sum_ris_weights += ris_samples[1].ris_weight;
		num_ris_candidates++;
	}
	assert(ris_samples[1].ris_weight >= 0.0f);
	assert(sum_ris_weights >= 0.0f);

	// ------------------------------------------------------------------------------
	// Sample/Select a sample from the RIS candidates proportional to the target
	// ------------------------------------------------------------------------------
	if (num_ris_candidates == 0 || !(sum_ris_weights > 1e-10f)) {
		*bsdf_pdf = 0.0f;
		*mis_pdf = 0.0f;
		return false;
	}

	float rand_ris_select = rand_bsdf_guiding * sum_ris_weights;

	float sum_ris = 0.0f;
	for (int i = 0; i < 2; i++) {
		sum_ris += ris_samples[i].ris_weight;
		if (rand_ris_select <= sum_ris) {
			ris_idx = i;
			break;
		}
	}

	assert(sum_ris >= 0.0f);
	assert(ris_idx < 2);

	// ------------------------------------------------------------------------------
	// Fill in the sample data for the selected RIS candidate
	// ------------------------------------------------------------------------------

	//ris_idx = 0;

	guide_pdf = ris_samples[ris_idx].ris_target * (2.0f / sum_ris_weights);
	*unguided_bsdf_pdf = ris_samples[ris_idx].bsdf_pdf;
	*mis_pdf = 0.5f * (ris_samples[ris_idx].bsdf_pdf + ris_samples[ris_idx].guide_pdf);
	*bsdf_pdf = guide_pdf;

	//*wo = ris_samples[ris_idx].wo;
	//label = ris_samples[ris_idx].label;

	*sampled_roughness = ris_samples[ris_idx].sampled_roughness;
	*eta = ris_samples[ris_idx].eta;
	*bsdf_eval = ris_samples[ris_idx].bsdf_eval;
	
	//localSampledDir = Vector(ris_samples[ris_idx].wo.x, ris_samples[ris_idx].wo.y, ris_samples[ris_idx].wo.z);
	//
	////*absCosSampledDir = fabsf(CosTheta(localSampledDir));
	//// to world
	//Vector ww = bsdf->GetFrame().ToWorld(localSampledDir);
	//ww = Normalize(ww);
	//*wo = float3(ww.x, ww.y, ww.z);
	*wo = ris_samples[ris_idx].wo;

	//assert(isfinite_safe(guide_pdf));
	//assert(isfinite_safe(*bsdf_pdf));

	if (!(*bsdf_pdf > 1e-10f)) {
		*bsdf_pdf = 0.0f;
		*mis_pdf = 0.0f;
		return  false;
	}

	assert(*bsdf_pdf > 0.0f);
	assert(*bsdf_pdf >= 1e-20f);
	assert(guide_pdf >= 0.0f);

	/// select label sampled_roughness and eta
	//if (ris_idx == 1 && ris_samples[1].bsdf_pdf > 0.0f) {
	//	float rnd = pathGuiding->path_state_rng_1D();
	//	//float rnd = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_RIS_GUIDING_1);

	//	float sum_pdfs = 0.0f;
	//	int idx = -1;
	//	for (int i = 0; i < sd.num_closure; i++) {
	//		sum_pdfs += unguided_bsdf_pdfs[i];
	//		if (rnd <= sum_pdfs) {
	//			idx = i;
	//			break;
	//		}
	//	}
	//	// assert(idx >= 0);
	//	/* Set the default idx to the last in the list.
	//	 * in case of numerical problems and rand_bsdf_guiding is just >=1.0f and
	//	 * the sum of all unguided_bsdf_pdfs is just < 1.0f. */
	//	idx = (rnd > sum_pdfs) ? sd.num_closure - 1 : idx;

	//	//label =  bsdf_label(kg, &sd->closure[idx], *wo);
	//	///bsdf_roughness_eta(kg, &sd->closure[idx], sampled_roughness, eta);
	//}

	/*assert(isfinite_safe(*bsdf_pdf));
	assert(*bsdf_pdf >= 0.0f);
	assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);*/

	if (!(*unguided_bsdf_pdf > 0.0f)) {
		*bsdf_pdf = 0.0f;
		*mis_pdf = 0.0f;
	}

	return true;
}

// ========================================================

bool PathTracer::SampleBSDF(BSDF* bsdf, Vector* sampledDir,
	const float u0, const float u1,
	Spectrum* bsdf_weight, // bsdf_eval / bsdf_pdf
	float* bsdf_pdf, 
	float* absCosSampledDir,
	BSDFEvent* event, PathGuiding* pathGuiding) const
{
	const bool use_surface_guiding = pathGuiding->state->use_surface_guiding;
	const float guiding_sampling_prob = pathGuiding->state->surface_guiding_sampling_prob;

	//   if !(use_surface_guiding && guiding_sampling_prob > 0.0f) {
	if (!use_surface_guiding || guiding_sampling_prob == 0.0f)
	{
		HitPoint& hitPoint = bsdf->hitPoint;

		Vector localFixedDir = bsdf->GetFrame().ToLocal(hitPoint.fixedDir);
		Vector localSampledDir;

		Spectrum result = bsdf->GetMaterial()->Sample(hitPoint,
			localFixedDir, &localSampledDir, u0, u1, hitPoint.passThroughEvent,
			bsdf_pdf, event);
		if (result.Black())
			return false;

		// =============================

		const Vector& localLightDir = hitPoint.fromLight ? localFixedDir : localSampledDir;
		const Vector& localEyeDir = hitPoint.fromLight ? localSampledDir : localFixedDir;

		float directPDF;
		Spectrum res1 = bsdf->GetMaterial()->Evaluate(hitPoint, localLightDir, localEyeDir, event,
			&directPDF);

		res1 = res1 / directPDF;
		int u = 0;

	/*	if ((res1 - result).Abs().Max() > 0.01)
		{
			int y = 0;
		}
		if (abs(directPDF - *bsdf_pdf) > 0.01)
		{
			int y = 0;
		}*/
		// 
		// =============

		*absCosSampledDir = fabsf(CosTheta(localSampledDir));
		*sampledDir = bsdf->GetFrame().ToWorld(localSampledDir);

		// Shadow terminator artefact avoidance
		if ((*event & REFLECT) &&
			(*event & (DIFFUSE | GLOSSY))
			&& (hitPoint.shadeN != hitPoint.interpolatedN)) {
			const Vector& lightDir = hitPoint.fromLight ? hitPoint.fixedDir : (*sampledDir);

			result *= _ShadowTerminatorAvoidanceFactor(hitPoint.GetLandingInterpolatedN(),
				hitPoint.GetLandingShadeN(), lightDir);
		}

		// Adjoint BSDF
		if (hitPoint.fromLight) {
			const float absDotFixedDirNG = AbsDot(hitPoint.fixedDir, hitPoint.geometryN);
			const float absDotSampledDirNG = AbsDot(*sampledDir, hitPoint.geometryN);
			result *= (absDotSampledDirNG / absDotFixedDirNG);
		}
		*bsdf_weight = result;
		return true;
	}
	else
	{
		BsdfEval bsdf_eval;
		float3 wo;
		//float bsdf_pdf;
		float mis_pdf;
		float unguided_bsdf_pdf;
		float2 sampled_roughness;
		float eta;

		bool ok = SampleBSDF_guiding(pathGuiding,bsdf, u0, u1,  event,
			&bsdf_eval, &wo, bsdf_pdf, &mis_pdf, &unguided_bsdf_pdf, &sampled_roughness, &eta);
		if (ok)
		{
			

			HitPoint& hitPoint = bsdf->hitPoint;
			*bsdf_weight = bsdf_eval.sum / *bsdf_pdf;
			*sampledDir = Vector(wo.x, wo.y, wo.z);

			auto localSampledDir = hitPoint.GetFrame().ToLocal(*sampledDir);
			*absCosSampledDir = fabsf(CosTheta(localSampledDir));
			//*pdfW = bsdf_pdf;

			// Shadow terminator artefact avoidance
			if ((*event & REFLECT) &&
				(*event & (DIFFUSE | GLOSSY))
				&& (hitPoint.shadeN != hitPoint.interpolatedN)) {
				const Vector& lightDir = hitPoint.fromLight ? hitPoint.fixedDir : (*sampledDir);

				(*bsdf_weight) *= _ShadowTerminatorAvoidanceFactor(hitPoint.GetLandingInterpolatedN(),
					hitPoint.GetLandingShadeN(), lightDir);
			}

			// Adjoint BSDF
			if (hitPoint.fromLight) {
				const float absDotFixedDirNG = AbsDot(hitPoint.fixedDir, hitPoint.geometryN);
				const float absDotSampledDirNG = AbsDot(*sampledDir, hitPoint.geometryN);
				(*bsdf_weight) *= (absDotSampledDirNG / absDotFixedDirNG);
			}
			return true;
		}
		else
			return false;
	}
}