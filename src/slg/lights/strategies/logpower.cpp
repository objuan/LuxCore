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

#include "slg/lights/strategies/logpower.h"
#include "slg/scene/scene.h"

using namespace std;
using namespace luxrays;
using namespace slg;

//------------------------------------------------------------------------------
// LightStrategyLogPower
//------------------------------------------------------------------------------

void LightStrategyLogPower::Preprocess(const Scene *scn, const LightStrategyTask taskType,
			const bool useRTMode) {
	// Delete old lightsDistribution
	delete lightsDistribution[0];
	delete lightsDistribution[1];
	delete lightsDistribution[2];
	lightsDistribution[0] = nullptr;
	lightsDistribution[1] = nullptr;
	lightsDistribution[2] = nullptr;

	DistributionLightStrategy::Preprocess(scn, taskType);

	const u_int lightCount = scene->lightDefs.GetSize();
	if (lightCount == 0)
		return;

	vector<float> lightPower_GROUP_0;
	vector<float> lightPower_GROUP_1;
	vector<float> lightPower_ALL;

	lightPower_GROUP_0.reserve(lightCount);
	lightPower_GROUP_1.reserve(lightCount);
	lightPower_ALL.reserve(lightCount);

	const vector<LightSource *> &lights = scene->lightDefs.GetLightSources();
	for (u_int i = 0; i < lightCount; ++i) {
		const LightSource *l = lights[i];
		float power = logf(1.f + l->GetPower(*scene)) * l->GetImportance();
		int lightGroup = l->GetID()==0 ? 0 : 1;
		bool g0 = lightGroup == 0; //l->IsIntersectable();

		switch (taskType) {
			case TASK_EMIT: {
				lightPower_ALL.push_back(power);

				if (g0) lightPower_GROUP_0.push_back(power); else lightPower_GROUP_0.push_back(0.f);
				if (!g0) lightPower_GROUP_1.push_back(power); else lightPower_GROUP_1.push_back(0.f);

				break;
			}
			case TASK_ILLUMINATE: {
				if (l->IsDirectLightSamplingEnabled()){
					lightPower_ALL.push_back(power);
					if (g0)  lightPower_GROUP_0.push_back(power); else lightPower_GROUP_0.push_back(0.f);
					if (!g0)  lightPower_GROUP_1.push_back(power); else lightPower_GROUP_1.push_back(0.f);
				}
				else
				{
					lightPower_ALL.push_back(0.f);
					lightPower_GROUP_0.push_back(0.f);
					lightPower_GROUP_1.push_back(0.f);
				}
				
				break;
			}
			case TASK_INFINITE_ONLY: {
				if (l->IsInfinite())
				{
					lightPower_ALL.push_back(power);
					if (g0) lightPower_GROUP_0.push_back(power); else lightPower_GROUP_0.push_back(0.f);
					if (!g0) lightPower_GROUP_1.push_back(power); else lightPower_GROUP_1.push_back(0.f);
				}
				else
				{
					lightPower_ALL.push_back(0.f);
					lightPower_GROUP_0.push_back(0.f);
					lightPower_GROUP_1.push_back(0.f);
				}
				break;
			}
			default:
				throw runtime_error("Unknown task in LightStrategyLogPower::Preprocess(): " + ToString(taskType));
		}
	}

	// Build the data to power based light sampling
	lightsDistribution[0] = new Distribution1D(&lightPower_ALL[0], lightCount);
	lightsDistribution[1] = new Distribution1D(&lightPower_GROUP_0[0], lightCount);
	lightsDistribution[2] = new Distribution1D(&lightPower_GROUP_1[0], lightCount);
}

// Static methods used by LightStrategyRegistry

Properties LightStrategyLogPower::ToProperties(const Properties &cfg) {
	return Properties() <<
			cfg.Get(GetDefaultProps().Get("lightstrategy.type"));
}

LightStrategy *LightStrategyLogPower::FromProperties(const Properties &cfg) {
	return new LightStrategyLogPower();
}

const Properties &LightStrategyLogPower::GetDefaultProps() {
	static Properties props = Properties() <<
			LightStrategy::GetDefaultProps() <<
			Property("lightstrategy.type")(GetObjectTag());

	return props;
}
