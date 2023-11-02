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

#ifndef LIBXL_ICONDFORMATT_H
#define LIBXL_ICONDFORMATT_H

#include "setup.h"
#include "enum.h"

namespace libxl
{
    template<class TCHAR> struct IConditionalFormatT;

    template<class TCHAR>
    struct IConditionalFormattingT
    {
        virtual void XLAPIENTRY addRange(int rowFirst, int rowLast, int colFirst, int colLast) = 0;
        virtual void XLAPIENTRY addRule(CFormatOperator op, IConditionalFormatT<TCHAR>* cFormat, const TCHAR* value1, const TCHAR* value2 = 0, CFormatType type = CFORMAT_CELLIS) = 0;
    };
}

#endif


