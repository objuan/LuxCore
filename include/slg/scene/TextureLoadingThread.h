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

#ifndef _SLG_TextureLoadingThread
#define	_SLG_TextureLoadingThread

#include <string>
#include <iostream>
#include <fstream>

#include "luxrays/core/intersectiondevice.h"
#include "luxrays/core/accelerator.h"
#include "luxrays/core/geometry/transform.h"
#include "luxrays/core/geometry/motionsystem.h"
#include "luxrays/utils/mc.h"
#include "luxrays/utils/mcdistribution.h"
#include "luxrays/utils/properties.h"
#include "luxrays/utils/serializationutils.h"
//#include "slg/core/sdl.h"
//#include "slg/cameras/camera.h"
//#include "slg/editaction.h"
//#include "slg/lights/light.h"
//#include "slg/lights/lightsourcedefs.h"
//#include "slg/shapes/strands.h"
//#include "slg/textures/texture.h"
//#include "slg/textures/texturedefs.h"
//#include "slg/textures/mapping/mapping.h"
//#include "slg/materials/materialdefs.h"
//#include "slg/bsdf/bsdf.h"
//#include "slg/volumes/volume.h"
//#include "slg/scene/sceneobjectdefs.h"
//#include "slg/scene/extmeshcache.h"
//#include "slg/scene/colorspaceconverters.h"

namespace slg {

	class Scene;

	class TextureLoadingTask
	{
		int id;
		slg::Scene* scene;
		//luxrays::Properties props;
		std::string props;

	public:
		std::string txtName;

		TextureLoadingTask(int id, slg::Scene* scene,  const luxrays::Properties& props, const std::string& txtName) :id(id), scene(scene), txtName(txtName) {
			this->props = props.ToString();
		}

		void load(int threadId);
		
	};

	// ===========================

	class TextureLoadingThread
	{
		#define threadCount 8
		boost::mutex mtx_;
		boost::thread* thread[threadCount];
		boost::mutex txtPathsMutex;
		std::string txtPaths[threadCount];
		slg::Scene* scene;
		std::vector< TextureLoadingTask *> taskList;
		bool ancora;
		int waitingCount;

		double lastUpdate;
	public:
		TextureLoadingThread(slg::Scene* scene) ;
		~TextureLoadingThread();

		void push(const luxrays::Properties& props, const std::string& txtName);

		void run(int threadId);

		void process();

		void clear() {
			taskList.clear();
		}

		void waitEnd();

		void lockTexturePath(const std::string &txtPath, int threadId);
		void unlockTexturePath(int threadId);

		std::vector<slg::Texture*>* toUpdate;
	};
}
#endif	/* _SLG_SCENE_H */
