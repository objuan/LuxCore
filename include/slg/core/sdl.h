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

#ifndef _SLG_SDL_H
#define	_SLG_SDL_H

#include <sstream>

#include "luxrays/luxrays.h"
#include "luxrays/core/color/color.h"

namespace slg {

extern void (*SLG_SDLDebugHandler)(const char *msg);

#define SDL_LOG(a) { if (slg::SLG_SDLDebugHandler) { std::stringstream _LR_LOG_LOCAL_SS; time_t now = time(NULL); struct tm t = *localtime(&now); _LR_LOG_LOCAL_SS << "[" << t.tm_mday << "-" << t.tm_mon + 1 << "-" << t.tm_year + 1900 << " " << t.tm_hour << ":" << t.tm_min << ":" << t.tm_sec << "] ";  _LR_LOG_LOCAL_SS << a; slg::SLG_SDLDebugHandler(_LR_LOG_LOCAL_SS.str().c_str()); } }

/*std::string GetTimeStamp_log()
{
	time_t now = time(NULL);
	struct tm t = *localtime(&now);
	std::stringstream ss;
	ss << t.tm_mday << "-" << t.tm_mon + 1 << "-" << t.tm_year + 1900 << " " << t.tm_hour << ":" << t.tm_min << ":" << t.tm_sec;
	return ss.str();
}*/
}

#endif	/* _SLG_SDL_H */
