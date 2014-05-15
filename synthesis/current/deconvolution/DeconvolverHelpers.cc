/// @file DeconvolverHelpers.cc
/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///

#include <askap_synthesis.h>

#include <askap/AskapLogging.h>
#include <askap/AskapError.h>
#include <casa/OS/Directory.h>
#include <images/Images/PagedImage.h>
#include <casa/BasicSL/String.h>   // for downcase

ASKAP_LOGGER(logger, ".deconvolver.helpers");

#include <deconvolution/DeconvolverHelpers.h>

namespace askap {
    namespace synthesis {

        Array<Float> DeconvolverHelpers::getArrayFromImage(const String name,
                                                           const LOFAR::ParameterSet& parset)
        {
            String imageName = parset.getString(name, "");
            casa::Array<Float> imageArray;
            if (imageName != "") {
                casa::PagedImage<float> im(imageName);
                im.get(imageArray, True);
                ASKAPLOG_INFO_STR(logger, "Read image " << imageName << " into array " << name);
            } else {
                ASKAPLOG_INFO_STR(logger, "Image " << name << " not specified in parset");
            }
            return imageArray;
        }

        void DeconvolverHelpers::putArrayToImage(const casa::Array<Float> imageArray, const String name,
                const String templateName, const LOFAR::ParameterSet &parset)
        {
            String templateFile = parset.getString(templateName, templateName);
            String imageFile = parset.getString(name, name);
            {
                ASKAPLOG_INFO_STR(logger, "Writing array " << name << " into image " << imageFile);
                casa::IPosition minPos;
                casa::IPosition maxPos;
                Float minVal, maxVal;
                casa::minMax(minVal, maxVal, minPos, maxPos, imageArray);
                ASKAPLOG_INFO_STR(logger, "Maximum =  " << maxVal << " at location " << maxPos);
                ASKAPLOG_INFO_STR(logger, "Minimum = " << minVal << " at location " << minPos);
                ASKAPLOG_INFO_STR(logger, "Sum     = " << sum(imageArray));
            }
            //
            {
                {
                    casa::Directory tFile(templateFile);
                    tFile.copy(imageFile);
                }
                casa::PagedImage<Float> im(imageFile);
                ASKAPLOG_INFO_STR(logger, "Array shape " << imageArray.shape());
                ASKAPLOG_INFO_STR(logger, "Image shape " << im.shape());

                im.put(imageArray);
                im.flush();
            }
        }
    }

}
