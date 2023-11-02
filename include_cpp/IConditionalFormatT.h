///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                    LibXL C++ headers version 4.0.3                        //
//                                                                           //
//                 Copyright (c) 2008 - 2022 XLware s.r.o.                   //
//                                                                           //
//   THIS FILE AND THE SOFTWARE CONTAINED HEREIN IS PROVIDED 'AS IS' AND     //
//                COMES WITH NO WARRANTIES OF ANY KIND.                      //
//                                                                           //
//          Please define LIBXL_STATIC variable for static linking.          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBXL_ICFORMATT_H
#define LIBXL_ICFORMATT_H

#include "setup.h"
#include "enum.h"

namespace libxl
{
    template<class TCHAR>
    struct IConditionalFormatT
    {
        virtual    FillPattern XLAPIENTRY fillPattern() const = 0;
        virtual           void XLAPIENTRY setFillPattern(FillPattern pattern) = 0;

        virtual          Color XLAPIENTRY patternForegroundColor() const = 0;
        virtual           void XLAPIENTRY setPatternForegroundColor(Color color) = 0;

        virtual          Color XLAPIENTRY patternBackgroundColor() const = 0;
        virtual           void XLAPIENTRY setPatternBackgroundColor(Color color) = 0;

        virtual                           ~IConditionalFormatT() {}
    };

}

#endif

