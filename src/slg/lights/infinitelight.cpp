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

#include <memory>

#include "slg/bsdf/bsdf.h"
#include "slg/scene/scene.h"
#include "slg/lights/infinitelight.h"

using namespace std;
using namespace luxrays;
using namespace slg;


//------------------------------------------------------------------------------
// InfiniteLight
//------------------------------------------------------------------------------

InfiniteLight::InfiniteLight() :
	imageMap(NULL), imageMapWork(NULL),imageMapDistribution(nullptr), visibilityMapCache(nullptr) {
}

InfiniteLight::~InfiniteLight() {
	delete imageMapDistribution;
	delete visibilityMapCache;
	delete imageMapWork;
}

float GetLum(Spectrum& color) {
	return 0.212671f * color.c[0] + 0.715160f * color.c[1] + 0.072169f * color.c[2];
	//return 0.11f * color.c[0] + 0.59f * color.c[1] + 0.3f * color.c[2]; //GRAY
}

void AddLum(Spectrum& color,float factor) {
	color.c[0] = max(0.f,color.c[0] + (factor ));
	color.c[1] = max(0.f, color.c[1] + (factor ));
	color.c[2] = max(0.f, color.c[2] + (factor ));

	/*color.c[0] = color.c[0] + (factor * 0.11f);
	color.c[1] = color.c[1] + (factor * 0.59f);
	color.c[2] = color.c[2] + (factor * 0.3f);*/
}

void InfiniteLight::Preprocess() {
	EnvLightSource::Preprocess();


	SDL_LOG("InfiniteLisht preprocess. contrast " << contrast << " hue " << hue << " saturation " << saturation);

	imageMapWork = imageMap->Copy();
		
	ImageMapStorage *imageMapStorage = imageMap->GetStorage();
	ImageMapStorage* imageMapStorageDest = imageMapWork->GetStorage();

	//	this->gain = Spectrum(contrast);
//	this->adder = Spectrum(abs(brightness));

	
	float lum;
	float min_lum = 999999;
	float max_lum = -9999999;
	vector<float> data(imageMap->GetWidth() * imageMap->GetHeight());
	//float maxVal = -INFINITY;
	//float minVal = INFINITY;
	for (u_int y = 0; y < imageMap->GetHeight(); ++y) {
		for (u_int x = 0; x < imageMap->GetWidth(); ++x) {
			const u_int index = x + y * imageMap->GetWidth();

			if (sampleUpperHemisphereOnly && (y > imageMap->GetHeight() / 2))
				data[index] = 0.f;
			else
				data[index] = imageMapStorage->GetFloat(index);

			if (!IsValid(data[index]))
				throw runtime_error("Pixel (" + ToString(x) + ", " + ToString(y) + ") in infinite light has an invalid value: " + ToString(data[index]));

			lum = imageMapStorage->GetSpectrum(index).Y();

			min_lum = min(min_lum, lum);
			max_lum = max(max_lum, lum);
			//maxVal = Max(data[index], maxVal);
			//minVal = Min(data[index], minVal);
		}
	}
	
	//this->adjImageMap = new luxrays::Spectrum[imageMap->GetWidth() * imageMap->GetHeight()];

	float lumAll = max_lum - min_lum;
	float cutMin = (inBlack/100.f); // 0-1
	float cutMax = (inWhite / 100.f);  // 0-1
	float cutSize = cutMax - cutMin;

	float outMin = (outBlack / 100.f); // 0-1
	float outMax = (outWhite / 100.f);  // 0-1
	float outSize = outMax - outMin;
	float cc = 1 + 1 - cutSize;

	float lum_in,lum_out;

	const float brightness = 0.0f;
	const float a = 1.f + contrast;
	const float b = brightness - contrast * 0.5f;

	for (u_int y = 0; y < imageMap->GetHeight(); ++y) {
		for (u_int x = 0; x < imageMap->GetWidth(); ++x) {
			const u_int index = x + y * imageMap->GetWidth();

			// luminosity
			Spectrum sp = imageMapStorage->GetSpectrum(index);

/*
				for (int i = 0; i < 3; i++)
				{
					lum = sp.c[i];
					// scale
					lum -= 0.5f;
					lum *= cc;
					lum += 0.5f;

					lum_in = min(99999.f, max(0.f, lum));
					sp.c[i] = lum_in;


				}
				*/

			sp = sp * a + Spectrum(b);

			if (saturation != 1.0f || hue != 0.5f)
			{
				Spectrum hsv = RgbToHsv(sp);
				if (hue != 0.5f)
				{
					hsv.c[0] += hue + .5f;
					hsv.c[0] = fmodf(hsv.c[0], 1.f);
				}
				if (saturation != 1.0f)
					hsv.c[1] *= saturation;
				sp = HsvToRgb(hsv);
			}

			for (int i = 0; i < 3; i++)
			{
				if (sp.c[i] < 0)
					sp.c[i] = 0;
			}


			imageMapStorageDest->SetSpectrum(index, sp);
		}
	}

	//SLG_LOG("InfiniteLight luminance  Max=" << maxVal << " Min=" << minVal);

	imageMapDistribution = new Distribution2D(&data[0], imageMap->GetWidth(), imageMap->GetHeight());
}

void InfiniteLight::GetPreprocessedData(const Distribution2D **imageMapDistributionData,
		const EnvLightVisibilityCache **elvc) const {
	if (imageMapDistributionData)
		*imageMapDistributionData = imageMapDistribution;
	if (elvc)
		*elvc = visibilityMapCache;
}

float InfiniteLight::GetPower(const Scene &scene) const {
	const float envRadius = GetEnvRadius(scene);

	// TODO: I should consider sampleUpperHemisphereOnly here
	return temperatureScale.Y() * gain.Y() * imageMap->GetSpectrumMeanY() *
			(4.f * M_PI * M_PI * envRadius * envRadius);
}

UV InfiniteLight::GetEnvUV(const luxrays::Vector &dir) const {
	UV uv;
	const Vector localDir = Normalize(Inverse(lightToWorld) * -dir);
	ToLatLongMapping(localDir, &uv.u, &uv.v);
	
	return uv;
}

Spectrum InfiniteLight::GetRadiance(const Scene &scene,
		const BSDF *bsdf, const Vector &dir,
		float *directPdfA, float *emissionPdfW) const {
	const Vector localDir = Normalize(Inverse(lightToWorld) * -dir);

	float u, v, latLongMappingPdf;
	ToLatLongMapping(localDir, &u, &v, &latLongMappingPdf);
	if (latLongMappingPdf == 0.f)
		return Spectrum();

	const float distPdf = imageMapDistribution->Pdf(u, v);
	if (directPdfA) {
		if (!bsdf)
			*directPdfA = 0.f;
		else if (visibilityMapCache && visibilityMapCache->IsCacheEnabled(*bsdf)) {
			*directPdfA = visibilityMapCache->Pdf(*bsdf, u, v) * latLongMappingPdf;
		} else
			*directPdfA = distPdf * latLongMappingPdf;
	}

	if (emissionPdfW) {
		const float envRadius = GetEnvRadius(scene);
		*emissionPdfW = distPdf * latLongMappingPdf / (M_PI * envRadius * envRadius);
	}

	Spectrum result =  temperatureScale * gain * imageMapWork->GetSpectrum(UV(u, v)) ;
	
	//else if (brightness < 0) result = result - adder;
	return result;
}

Spectrum InfiniteLight::Emit(const Scene &scene,
		const float time, const float u0, const float u1,
		const float u2, const float u3, const float passThroughEvent,
		Ray &ray, float &emissionPdfW,
		float *directPdfA, float *cosThetaAtLight) const {
	float uv[2];
	float distPdf;
	imageMapDistribution->SampleContinuous(u0, u1, uv, &distPdf);
	if (distPdf == 0.f)
		return Spectrum();
	
	Vector localDir;
	float latLongMappingPdf;
	FromLatLongMapping(uv[0], uv[1], &localDir, &latLongMappingPdf);
	if (latLongMappingPdf == 0.f)
		return Spectrum();

	// Compute the ray direction
	const Vector rayDir = -Normalize(lightToWorld * localDir);

	// Compute the ray origin
	Vector x, y;
    CoordinateSystem(-rayDir, &x, &y);
    float d1, d2;
    ConcentricSampleDisk(u2, u3, &d1, &d2);

	const Point worldCenter = scene.dataSet->GetBSphere().center;
	const float envRadius = GetEnvRadius(scene);
	const Point pDisk = worldCenter + envRadius * (d1 * x + d2 * y);
	const Point rayOrig = pDisk - envRadius * rayDir;

	// Compute InfiniteLight ray weight
	emissionPdfW = distPdf * latLongMappingPdf / (M_PI * envRadius * envRadius);

	if (directPdfA)
		*directPdfA = distPdf * latLongMappingPdf;

	if (cosThetaAtLight)
		*cosThetaAtLight = Dot(Normalize(worldCenter - rayOrig), rayDir);

	Spectrum result = temperatureScale * gain * imageMap->GetSpectrum(uv);
	
	//else if (brightness < 0) result = result - adder;

	assert (!result.IsNaN() && !result.IsInf() && !result.IsNeg());

	ray.Update(rayOrig, rayDir, time);

	return result;
}

Spectrum InfiniteLight::Illuminate(const Scene &scene, const BSDF &bsdf,
		const float time, const float u0, const float u1, const float passThroughEvent,
        Ray &shadowRay, float &directPdfW,
		float *emissionPdfW, float *cosThetaAtLight) const {
	float uv[2];
	float distPdf;	
	if (visibilityMapCache && visibilityMapCache->IsCacheEnabled(bsdf))
		visibilityMapCache->Sample(bsdf, u0, u1, uv, &distPdf);
	else
		imageMapDistribution->SampleContinuous(u0, u1, uv, &distPdf);
	if (distPdf == 0.f)
		return Spectrum();

	Vector localDir;
	float latLongMappingPdf;
	FromLatLongMapping(uv[0], uv[1], &localDir, &latLongMappingPdf);
	if (latLongMappingPdf == 0.f)
		return Spectrum();

	const Vector shadowRayDir = Normalize(lightToWorld * localDir);
	
	const Point worldCenter = scene.dataSet->GetBSphere().center;
	const float envRadius = GetEnvRadius(scene);

	const Point shadowRayOrig = bsdf.GetRayOrigin(shadowRayDir);
	const Vector toCenter(worldCenter - shadowRayOrig);
	const float centerDistanceSquared = Dot(toCenter, toCenter);
	const float approach = Dot(toCenter, shadowRayDir);
	const float shadowRayDistance = approach + sqrtf(Max(0.f, envRadius * envRadius -
		centerDistanceSquared + approach * approach));

	const Point emisPoint(shadowRayOrig + shadowRayDistance * shadowRayDir);
	const Normal emisNormal(Normalize(worldCenter - emisPoint));

	const float cosAtLight = Dot(emisNormal, -shadowRayDir);
	if (cosAtLight < DEFAULT_COS_EPSILON_STATIC)
		return Spectrum();
	if (cosThetaAtLight)
		*cosThetaAtLight = cosAtLight;

	directPdfW = distPdf * latLongMappingPdf;
	assert (!isnan(directPdfW) && !isinf(directPdfW) && (directPdfW > 0.f));

	if (emissionPdfW)
		*emissionPdfW = distPdf * latLongMappingPdf / (M_PI * envRadius * envRadius);

	 Spectrum result = temperatureScale * gain * imageMap->GetSpectrum(UV(uv[0], uv[1])) ;
	
	// else if (brightness < 0) result = result - adder;
	 
	 assert (!result.IsNaN() && !result.IsInf() && !result.IsNeg());

	shadowRay = Ray(shadowRayOrig, shadowRayDir, 0.f, shadowRayDistance, time);

	return result;
}

void InfiniteLight::UpdateVisibilityMap(const Scene *scene, const bool useRTMode) {
	delete visibilityMapCache;
	visibilityMapCache = nullptr;

	if (useRTMode)
		return;

	if (useVisibilityMapCache) {
		// Scale the infinitelight image map to the requested size
		unique_ptr<ImageMap> luminanceMapImage(imageMap->Copy());
		// Select the image luminance
		luminanceMapImage->SelectChannel(ImageMapStorage::WEIGHTED_MEAN);
		luminanceMapImage->Preprocess();

		visibilityMapCache = new EnvLightVisibilityCache(scene, this,
				luminanceMapImage.get(), visibilityMapCacheParams);		
		visibilityMapCache->Build();
	}
}

Spectrum InfiniteLight::RgbToHsv(const Spectrum& rgb) {
	float cmax, cmin, h, s, v, cdelta;

	cmax = Max(rgb.c[0], Max(rgb.c[1], rgb.c[2]));
	cmin = Min(rgb.c[0], Min(rgb.c[1], rgb.c[2]));
	cdelta = cmax - cmin;

	v = cmax;

	if (cmax != 0.f)
		s = cdelta / cmax;
	else {
		s = 0.f;
		h = 0.f;
	}

	if (s != 0.0f) {
		Spectrum c;
		float icdelta = 1.f / cdelta;
		c.c[0] = (cmax - rgb.c[0]) * icdelta;
		c.c[1] = (cmax - rgb.c[1]) * icdelta;
		c.c[2] = (cmax - rgb.c[2]) * icdelta;

		if (rgb.c[0] == cmax)
			h = c.c[2] - c.c[1];
		else if (rgb.c[1] == cmax)
			h = 2.f + c.c[0] - c.c[2];
		else
			h = 4.f + c.c[1] - c.c[0];

		h /= 6.f;

		if (h < 0.f)
			h += 1.f;
	}
	else
		h = 0.f;

	return Spectrum(h, s, v);
}

Spectrum InfiniteLight::HsvToRgb(const Spectrum& hsv) {
	float i, f, p, q, t, h, s, v;

	h = hsv.c[0];
	s = hsv.c[1];
	v = hsv.c[2];

	if (s != 0.f) {
		if (h == 1.f)
			h = 0.f;

		h *= 6.f;
		i = Floor2Int(h);
		f = h - i;

		p = v * (1.f - s);
		q = v * (1.f - (s * f));
		t = v * (1.f - (s * (1.f - f)));

		if (i == 0.f) return Spectrum(v, t, p);
		else if (i == 1.f) return Spectrum(q, v, p);
		else if (i == 2.f) return Spectrum(p, v, t);
		else if (i == 3.f) return Spectrum(p, q, v);
		else if (i == 4.f) return Spectrum(t, p, v);
		else return Spectrum(v, p, q);
	}
	else
		return Spectrum(v, v, v);
}

Properties InfiniteLight::ToProperties(const ImageMapCache &imgMapCache, const bool useRealFileName) const {
	const string prefix = "scene.lights." + GetName();
	Properties props = EnvLightSource::ToProperties(imgMapCache, useRealFileName);

	props.Set(Property(prefix + ".type")("infinite"));
	props.Set(Property(prefix + ".hue")(hue));
	props.Set(Property(prefix + ".contrast")(contrast));
	props.Set(Property(prefix + ".saturation")(saturation));
	props.Set(Property(prefix + ".levelMinPerc")(0.f));
	props.Set(Property(prefix + ".levelMaxPerc")(100.f));
	props.Set(Property(prefix + ".gammaCorrection")(1.f));
	const string fileName = useRealFileName ?
		imageMap->GetName() : imgMapCache.GetSequenceFileName(imageMap);
	props.Set(Property(prefix + ".file")(fileName));
	props.Set(imageMap->ToProperties(prefix, false));
	props.Set(Property(prefix + ".gamma")(1.f));
	props.Set(Property(prefix + ".sampleupperhemisphereonly")(sampleUpperHemisphereOnly));

	props.Set(Property(prefix + ".visibilitymapcache.enable")(useVisibilityMapCache));
	if (useVisibilityMapCache)
		props << EnvLightVisibilityCache::Params2Props(prefix, visibilityMapCacheParams);

	return props;
}
