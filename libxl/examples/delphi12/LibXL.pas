///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                  LibXL Delphi wrapper version 3.9.4                       //
//                                                                           //
//       Copyright (c) 2008 - 2021 Dmytro Skrypnyk and XLware s.r.o.         //
//                         All rights reserved.                              //
//                                                                           //
// This file forms part of LibXL product. It may not be used independenty.   //
// See the LibXL licence agreement for further restrictions.                 //
//                                                                           //
// THIS FILE AND THE SOFTWARE CONTAINED HEREIN IS PROVIDED 'AS IS' AND       //
// COMES WITH NO WARRANTIES OF ANY KIND.                                     //
//                                                                           //
// XLWARE SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE             //
// OR ANYTHING DUE TO THE USE OF THIS FILE HOWEVER THEY MAY BE CAUSED.       //
//                                                                           //
// For more information on the LibXL, go to: http://www.libxl.com            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

unit LibXL;

{$Z4}

interface

uses SysUtils;

type
  AlignH = (ALIGNH_GENERAL, ALIGNH_LEFT, ALIGNH_CENTER, ALIGNH_RIGHT, ALIGNH_FILL, ALIGNH_JUSTIFY,
            ALIGNH_MERGE, ALIGNH_DISTRIBUTED);
type
  AlignV = (ALIGNV_TOP, ALIGNV_CENTER, ALIGNV_BOTTOM, ALIGNV_JUSTIFY, ALIGNV_DISTRIBUTED);

type
  BorderStyle = (BORDERSTYLE_NONE, BORDERSTYLE_THIN, BORDERSTYLE_MEDIUM, BORDERSTYLE_DASHED,
                 BORDERSTYLE_DOTTED, BORDERSTYLE_THICK, BORDERSTYLE_DOUBLE, BORDERSTYLE_HAIR,
                 BORDERSTYLE_MEDIUMDASHED, BORDERSTYLE_DASHDOT, BORDERSTYLE_MEDIUMDASHDOT,
                 BORDERSTYLE_DASHDOTDOT, BORDERSTYLE_MEDIUMDASHDOTDOT, BORDERSTYLE_SLANTDASHDOT);
type
  BorderDiagonal = (BORDERDIAGONAL_NONE, BORDERDIAGONAL_DOWN, BORDERDIAGONAL_UP,
                    BORDERDIAGONAL_BOTH);
type
  FillPattern = (FILLPATTERN_NONE, FILLPATTERN_SOLID, FILLPATTERN_GRAY50, FILLPATTERN_GRAY75,
                 FILLPATTERN_GRAY25, FILLPATTEN_HORSTRIPE, FILLPATTERN_VERSTRIPE,
                 FILLPATTERN_REVDIAGSTRIPE, FILLPATTERN_DIAGSTRIPE, FILLPATTERN_DIAGCROSSHATCH,
                 FILLPATTERN_THICKDIAGCROSSHATCH, FILLPATTERN_THINHORSTRIPE,
                 FILLPATTERN_THINVERSTRIPE, FILLPATTERN_THINREVDIAGSTRIPE,
                 FILLPATTERN_THINDIAGSTRIPE, FILLPATTERN_THINHORCROSSHATCH,
                 FILLPATTERN_THINDIAGCROSSHATCH, FILLPATTERN_GRAY12P5, FILLPATTERN_GRAY6P25);
type
  NumFormat = (NUMFORMAT_GENERAL, NUMFORMAT_NUMBER, NUMFORMAT_NUMBER_D2, NUMFORMAT_NUMBER_SEP,
               NUMFORMAT_NUMBER_SEP_D2, NUMFORMAT_CURRENCY_NEGBRA, NUMFORMAT_CURRENCY_NEGBRARED,
               NUMFORMAT_CURRENCY_D2_NEGBRA, NUMFORMAT_CURRENCY_D2_NEGBRARED, NUMFORMAT_PERCENT,
               NUMFORMAT_PERCENT_D2, NUMFORMAT_SCIENTIFIC_D2, NUMFORMAT_FRACTION_ONEDIG,
               NUMFORMAT_FRACTION_TWODIG, NUMFORMAT_DATE, NUMFORMAT_CUSTOM_D_MON_YY,
               NUMFORMAT_CUSTOM_D_MON, NUMFORMAT_CUSTOM_MON_YY, NUMFORMAT_CUSTOM_HMM_AM,
               NUMFORMAT_CUSTOM_HMMSS_AM, NUMFORMAT_CUSTOM_HMM, NUMFORMAT_CUSTOM_HMMSS,
               NUMFORMAT_CUSTOM_MDYYYY_HMM, NUMFORMAT_NUMBER_SEP_NEGBRA = 37,
               NUMFORMAT_NUMBER_SEP_NEGBRARED, NUMFORMAT_NUMBER_D2_SEP_NEGBRA,
               NUMFORMAT_NUMBER_D2_SEP_NEGBRARED, NUMFORMAT_ACCOUNT, NUMFORMAT_ACCOUNTCUR,
               NUMFORMAT_ACCOUNT_D2, NUMFORMAT_ACCOUNT_D2_CUR, NUMFORMAT_CUSTOM_MMSS,
               NUMFORMAT_CUSTOM_H0MMSS, NUMFORMAT_CUSTOM_MMSS0, NUMFORMAT_CUSTOM_000P0E_PLUS0,
               NUMFORMAT_TEXT);
type
  Color = (COLOR_BLACK = 8, COLOR_WHITE, COLOR_RED, COLOR_BRIGHTGREEN, COLOR_BLUE, COLOR_YELLOW,
           COLOR_PINK, COLOR_TURQUOISE, COLOR_DARKRED, COLOR_GREEN, COLOR_DARKBLUE,
           COLOR_DARKYELLOW, COLOR_VIOLET, COLOR_TEAL, COLOR_GRAY25, COLOR_GRAY50,
           COLOR_PERIWINKLE_CF, COLOR_PLUM_CF, COLOR_IVORY_CF, COLOR_LIGHTTURQUOISE_CF,
           COLOR_DARKPURPLE_CF, COLOR_CORAL_CF, COLOR_OCEANBLUE_CF, COLOR_ICEBLUE_CF,
           COLOR_DARKBLUE_CL, COLOR_PINK_CL, COLOR_YELLOW_CL, COLOR_TURQUOISE_CL, COLOR_VIOLET_CL,
           COLOR_DARKRED_CL, COLOR_TEAL_CL, COLOR_BLUE_CL, COLOR_SKYBLUE, COLOR_LIGHTTURQUOISE,
           COLOR_LIGHTGREEN, COLOR_LIGHTYELLOW, COLOR_PALEBLUE, COLOR_ROSE, COLOR_LAVENDER,
           COLOR_TAN, COLOR_LIGHTBLUE, COLOR_AQUA, COLOR_LIME, COLOR_GOLD, COLOR_LIGHTORANGE,
           COLOR_ORANGE, COLOR_BLUEGRAY, COLOR_GRAY40, COLOR_DARKTEAL, COLOR_SEAGREEN,
           COLOR_DARKGREEN, COLOR_OLIVEGREEN, COLOR_BROWN, COLOR_PLUM, COLOR_INDIGO, COLOR_GRAY80,
           COLOR_DEFAULT_FOREGROUND = $0040, COLOR_DEFAULT_BACKGROUND = $0041,
           COLOR_TOOLTIP = $0051, COLOR_AUTO = $7FFF);

type
  SheetType = (SHEETTYPE_SHEET, SHEETTYPE_CHART, SHEETTYPE_UNKNOWN);

type
  CellType = (CELLTYPE_EMPTY, CELLTYPE_NUMBER, CELLTYPE_STRING, CELLTYPE_BOOLEAN, CELLTYPE_BLANK,
              CELLTYPE_ERROR);
type
  ErrorType = (ERRORTYPE_NULL = $00, ERRORTYPE_DIV_0 = $07, ERRORTYPE_VALUE = $0F,
               ERRORTYPE_REF = $17, ERRORTYPE_NAME = $1D, ERRORTYPE_NUM = $24, ERRORTYPE_NA = $2A,
               ERRORTYPE_NOERROR = $FF);
type
  PictureType = (PICTURETYPE_PNG, PICTURETYPE_JPEG, PICTURETYPE_WMF, PICTURETYPE_DIB, PICTURETYPE_EMF,
                 PICTURETYPE_PICT, PICTURETYPE_TIFF, PICTURETYPE_ERROR = $FF);
type
  SheetState = (SHEETSTATE_VISIBLE, SHEETSTATE_HIDDEN, SHEETSTATE_VERYHIDDEN);
type
  Paper = (PAPER_DEFAULT, PAPER_LETTER, PAPER_LETTERSMALL, PAPER_TABLOID, PAPER_LEDGER, PAPER_LEGAL,
           PAPER_STATEMENT, PAPER_EXECUTIVE, PAPER_A3, PAPER_A4, PAPER_A4SMALL, PAPER_A5, PAPER_B4,
           PAPER_B5, PAPER_FOLIO, PAPER_QUATRO, PAPER_10x14, PAPER_10x17, PAPER_NOTE,
           PAPER_ENVELOPE_9, PAPER_ENVELOPE_10, PAPER_ENVELOPE_11, PAPER_ENVELOPE_12,
           PAPER_ENVELOPE_14, PAPER_C_SIZE, PAPER_D_SIZE, PAPER_E_SIZE, PAPER_ENVELOPE_DL,
           PAPER_ENVELOPE_C5, PAPER_ENVELOPE_C3, PAPER_ENVELOPE_C4, PAPER_ENVELOPE_C6,
           PAPER_ENVELOPE_C65, PAPER_ENVELOPE_B4, PAPER_ENVELOPE_B5, PAPER_ENVELOPE_B6,
           PAPER_ENVELOPE, PAPER_ENVELOPE_MONARCH, PAPER_US_ENVELOPE, PAPER_FANFOLD,
           PAPER_GERMAN_STD_FANFOLD, PAPER_GERMAN_LEGAL_FANFOLD, PAPER_B4_ISO,
           PAPER_JAPANESE_POSTCARD, PAPER_9x11, PAPER_10x11, PAPER_15x11, PAPER_ENVELOPE_INVITE,
           PAPER_US_LETTER_EXTRA = 50, PAPER_US_LEGAL_EXTRA, PAPER_US_TABLOID_EXTRA, PAPER_A4_EXTRA,
           PAPER_LETTER_TRANSVERSE, PAPER_A4_TRANSVERSE, PAPER_LETTER_EXTRA_TRANSVERSE,
           PAPER_SUPERA, PAPER_SUPERB, PAPER_US_LETTER_PLUS, PAPER_A4_PLUS, PAPER_A5_TRANSVERSE,
           PAPER_B5_TRANSVERSE, PAPER_A3_EXTRA, PAPER_A5_EXTRA, PAPER_B5_EXTRA, PAPER_A2,
           PAPER_A3_TRANSVERSE, PAPER_A3_EXTRA_TRANSVERSE, PAPER_JAPANESE_DOUBLE_POSTCARD, PAPER_A6,
           PAPER_JAPANESE_ENVELOPE_KAKU2, PAPER_JAPANESE_ENVELOPE_KAKU3,
           PAPER_JAPANESE_ENVELOPE_CHOU3, PAPER_JAPANESE_ENVELOPE_CHOU4, PAPER_LETTER_ROTATED,
           PAPER_A3_ROTATED, PAPER_A4_ROTATED, PAPER_A5_ROTATED, PAPER_B4_ROTATED, PAPER_B5_ROTATED,
           PAPER_JAPANESE_POSTCARD_ROTATED, PAPER_DOUBLE_JAPANESE_POSTCARD_ROTATED,
           PAPER_A6_ROTATED, PAPER_JAPANESE_ENVELOPE_KAKU2_ROTATED,
           PAPER_JAPANESE_ENVELOPE_KAKU3_ROTATED, PAPER_JAPANESE_ENVELOPE_CHOU3_ROTATED,
           PAPER_JAPANESE_ENVELOPE_CHOU4_ROTATED, PAPER_B6, PAPER_B6_ROTATED, PAPER_12x11,
           PAPER_JAPANESE_ENVELOPE_YOU4, PAPER_JAPANESE_ENVELOPE_YOU4_ROTATED, PAPER_PRC16K,
           PAPER_PRC32K, PAPER_PRC32K_BIG, PAPER_PRC_ENVELOPE1, PAPER_PRC_ENVELOPE2,
           PAPER_PRC_ENVELOPE3, PAPER_PRC_ENVELOPE4, PAPER_PRC_ENVELOPE5, PAPER_PRC_ENVELOPE6,
           PAPER_PRC_ENVELOPE7, PAPER_PRC_ENVELOPE8, PAPER_PRC_ENVELOPE9, PAPER_PRC_ENVELOPE10,
           PAPER_PRC16K_ROTATED, PAPER_PRC32K_ROTATED, PAPER_PRC32KBIG_ROTATED,
           PAPER_PRC_ENVELOPE1_ROTATED, PAPER_PRC_ENVELOPE2_ROTATED, PAPER_PRC_ENVELOPE3_ROTATED,
           PAPER_PRC_ENVELOPE4_ROTATED, PAPER_PRC_ENVELOPE5_ROTATED, PAPER_PRC_ENVELOPE6_ROTATED,
           PAPER_PRC_ENVELOPE7_ROTATED, PAPER_PRC_ENVELOPE8_ROTATED, PAPER_PRC_ENVELOPE9_ROTATED,
           PAPER_PRC_ENVELOPE10_ROTATED);
type
  Underline = (UNDERLINE_NONE, UNDERLINE_SINGLE, UNDERLINE_DOUBLE,
               UNDERLINE_SINGLEACC = $21, UNDERLINE_DOUBLEACC = $22);
type
  Script = (SCRIPT_NORMAL, SCRIPT_SUPER, SCRIPT_SUB);
type
  Scope = (SCOPE_UNDEFINED = -2, SCOPE_WORKBOOK = -1);
type
  Operator = (OPERATOR_EQUAL, OPERATOR_GREATER_THAN, OPERATOR_GREATER_THAN_OR_EQUAL, OPERATOR_LESS_THAN, OPERATOR_LESS_THAN_OR_EQUAL, OPERATOR_NOT_EQUAL );

type
  Filter = (FILTER_VALUE, FILTER_TOP10, FILTER_CUSTOM, FILTER_DYNAMIC, FILTER_COLOR, FILTER_ICON, FILTER_EXT, FILTER_NOT_SET);

type
  EnhancedProtection = (PROT_DEFAULT = -1, PROT_ALL = 0, PROT_OBJECTS = 1, PROT_SCENARIOS = 2, PROT_FORMAT_CELLS = 4, PROT_FORMAT_COLUMNS = 8, PROT_FORMAT_ROWS = 16,
                        PROT_INSERT_COLUMNS = 32, PROT_INSERT_ROWS = 64, PROT_INSERT_HYPERLINKS = 128, PROT_DELETE_COLUMNS = 256, PROT_DELETE_ROWS = 512,
                        PROT_SEL_LOCKED_CELLS = 1024, PROT_SORT = 2048, PROT_AUTOFILTER = 4096, PROT_PIVOTTABLES = 8192, PROT_SEL_UNLOCKED_CELLS = 16384);
type
  IgnoredError = (IERR_NO_ERROR = 0, IERR_EVAL_ERROR = 1, IERR_EMPTY_CELLREF = 2, IERR_NUMBER_STORED_AS_TEXT = 4, IERR_INCONSIST_RANGE = 8,
                  IERR_INCONSIST_FMLA = 16, IERR_TWODIG_TEXTYEAR = 32, IERR_UNLOCK_FMLA = 64, IERR_DATA_VALIDATION = 128);

type
  DataValidationType = (VALIDATION_TYPE_NONE, VALIDATION_TYPE_WHOLE, VALIDATION_TYPE_DECIMAL, VALIDATION_TYPE_LIST,
                        VALIDATION_TYPE_DATE, VALIDATION_TYPE_TIME, VALIDATION_TYPE_TEXTLENGTH, VALIDATION_TYPE_CUSTOM);

type
  DataValidationOperator = (VALIDATION_OP_BETWEEN, VALIDATION_OP_NOTBETWEEN, VALIDATION_OP_EQUAL, VALIDATION_OP_NOTEQUAL,
                            VALIDATION_OP_LESSTHAN, VALIDATION_OP_LESSTHANOREQUAL, VALIDATION_OP_GREATERTHAN, VALIDATION_OP_GREATERTHANOREQUAL);

type
  DataValidationErrorStyle = (VALIDATION_ERRSTYLE_STOP, VALIDATION_ERRSTYLE_WARNING, VALIDATION_ERRSTYLE_INFORMATION);

type
  CalcModeType = (CALCMODE_MANUAL, CALCMODE_AUTO, CALCMODE_AUTONOTABLE);

type
  BookHandle = Pointer;
type
  SheetHandle = Pointer;
type
  FormatHandle = Pointer;
type
  FontHandle = Pointer;
type
  RichStringHandle = Pointer;
type
  PFormatHandle = ^FormatHandle;
type
  PFontHandle = ^FontHandle;
type
  AutoFilterHandle = Pointer;
type
  FilterColumnHandle = Pointer;
type
  ByteArray = array of byte;
type
  PByteArray = ^ByteArray;
type
  PCardinal = ^Cardinal;
type
  PInteger = ^Integer;
type
  PBoolean = ^boolean;
type
  PDouble = ^double;

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlCreateBookC(): BookHandle cdecl;
external 'libxl' name 'xlCreateBookCW';

function xlCreateXMLBookC(): BookHandle cdecl;
external 'libxl' name 'xlCreateXMLBookCW';

function xlBookLoad(handle: BookHandle; filename: PWideChar): Integer cdecl;
external 'libxl' name 'xlBookLoadW';

function xlBookSave(handle: BookHandle; filename: PWideChar): Integer cdecl;
external 'libxl' name 'xlBookSaveW';

function xlBookLoadRaw(handle: BookHandle; data: PByteArray; size: Cardinal): Integer cdecl;
external 'libxl' name 'xlBookLoadRawW';

function xlBookSaveRaw(handle: BookHandle; data: PByteArray; size: PCardinal): Integer cdecl;
external 'libxl' name 'xlBookSaveRawW';

function xlBookAddSheet(handle: BookHandle; name: PWideChar; initSheet: SheetHandle): SheetHandle cdecl;
external 'libxl' name 'xlBookAddSheetW';

function xlBookInsertSheet(handle: BookHandle; index: Integer; name: PWideChar; initSheet: SheetHandle): SheetHandle cdecl;
external 'libxl' name 'xlBookInsertSheetW';

function xlBookGetSheet(handle: BookHandle; index: Integer): SheetHandle cdecl;
external 'libxl' name 'xlBookGetSheetW';

function xlBookSheetType(handle: BookHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlBookSheetTypeW';

function xlBookMoveSheet(handle: BookHandle; srcIndex: Integer; dstIndex: Integer): Integer cdecl;
external 'libxl' name 'xlBookMoveSheetW';

function xlBookDelSheet(handle: BookHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlBookDelSheetW';

function xlBookSheetCount(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookSheetCountW';

function xlBookAddFormat(handle: BookHandle; initFormat: FormatHandle): FormatHandle cdecl;
external 'libxl' name 'xlBookAddFormatW';

function xlBookAddFont(handle: BookHandle; initFont: FontHandle): FontHandle cdecl;
external 'libxl' name 'xlBookAddFontW';

function xlBookAddRichString(handle: BookHandle): RichStringHandle cdecl;
external 'libxl' name 'xlBookAddRichStringW';

function xlBookAddCustomNumFormat(handle: BookHandle; const customNumFormat: PWideChar): Integer cdecl;
external 'libxl' name 'xlBookAddCustomNumFormatW';

function xlBookCustomNumFormat(handle: BookHandle; fmt: Integer): PWideChar cdecl;
external 'libxl' name 'xlBookCustomNumFormatW';

function xlBookFormat(handle: BookHandle; index: Integer): FormatHandle cdecl;
external 'libxl' name 'xlBookFormatW';

function xlBookFormatSize(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookFormatSizeW';

function xlBookFont(handle: BookHandle; index: Integer): FontHandle cdecl;
external 'libxl' name 'xlBookFontW';

function xlBookFontSize(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookFontSizeW';

function xlBookDatePack(handle: BookHandle; year: Integer; month: Integer; day: Integer; hour: Integer; min: Integer; sec: Integer; msec: Integer): double cdecl;
external 'libxl' name 'xlBookDatePackW';

function xlBookDateUnpack(handle: BookHandle; value: double; year: PInteger; month: PInteger; day: PInteger; hour: PInteger; min: PInteger; sec: PInteger; msec: PInteger): Integer cdecl;
external 'libxl' name 'xlBookDateUnpackW';

function xlBookColorPack(handle: BookHandle; red: Integer; green: Integer; blue: Integer): Integer cdecl;
external 'libxl' name 'xlBookColorPackW';

procedure xlBookColorUnpack(handle: BookHandle; color: Integer; red: PInteger; green: PInteger; blue: PInteger) cdecl;
external 'libxl' name 'xlBookColorUnpackW';

function xlBookActiveSheet(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookActiveSheetW';

procedure xlBookSetActiveSheet(handle: BookHandle; index: Integer) cdecl;
external 'libxl' name 'xlBookSetActiveSheetW';

function xlBookPictureSize(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookPictureSizeW';

function xlBookGetPicture(hanlde: BookHandle; index: Integer; data: PByteArray; size: PCardinal): Integer cdecl;
external 'libxl' name 'xlBookGetPictureW';

function xlBookAddPicture(handle: BookHandle; filename: PWideChar): Integer cdecl;
external 'libxl' name 'xlBookAddPictureW';

function xlBookAddPicture2(handle: BookHandle; data: PByteArray; size: Cardinal): Integer cdecl;
external 'libxl' name 'xlBookAddPicture2W';

function xlBookAddPictureAsLink(handle: BookHandle; filename: PWideChar; insert: Integer): Integer cdecl;
external 'libxl' name 'xlBookAddPictureAsLinkW';

function xlBookDefaultFont(handle: BookHandle; fontSize: PInteger): PWideChar cdecl;
external 'libxl' name 'xlBookDefaultFontW';

procedure xlBookSetDefaultFont(handle: BookHandle; const name: PWideChar; fontSize: Integer) cdecl;
external 'libxl' name 'xlBookSetDefaultFontW';

function xlBookRefR1C1(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookRefR1C1W';

procedure xlBookSetRefR1C1(handle: BookHandle; refR1C1: Integer) cdecl;
external 'libxl' name 'xlBookSetRefR1C1W';

procedure xlBookSetKey(handle: BookHandle; const name: PWideChar; const key: PWideChar) cdecl;
external 'libxl' name 'xlBookSetKeyW';

function xlBookRgbMode(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookRgbModeW';

procedure xlBookSetRgbMode(handle: BookHandle; rgbMode: Integer) cdecl;
external 'libxl' name 'xlBookSetRgbModeW';

function xlBookCalcMode(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookCalcModeW';

procedure xlBookSetCalcMode(handle: BookHandle; calcMode: Integer) cdecl;
external 'libxl' name 'xlBookSetCalcModeW';

function xlBookIsDate1904(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookIsDate1904W';

procedure xlBookSetDate1904(handle: BookHandle; date1904: Integer) cdecl;
external 'libxl' name 'xlBookSetDate1904W';

function xlBookIsTemplate(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookIsTemplateW';

procedure xlBookSetTemplate(handle: BookHandle; tmpl: Integer) cdecl;
external 'libxl' name 'xlBookSetTemplateW';

function xlBookBiffVersion(handle: BookHandle): Integer cdecl;
external 'libxl' name 'xlBookBiffVersionW';

function xlBookSetLocale(handle: BookHandle; const locale: PAnsiChar): Integer cdecl;
external 'libxl' name 'xlBookSetLocaleW';

function xlBookErrorMessage(handle: BookHandle): PAnsiChar cdecl;
external 'libxl' name 'xlBookErrorMessageW';

procedure xlBookRelease(handle: BookHandle) cdecl;
external 'libxl' name 'xlBookReleaseW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlSheetCopyCell(handle: SheetHandle; rowSrc, colSrc, rowDst, colDst: Integer): Integer cdecl;
external 'libxl' name 'xlSheetCopyCellW';

function xlSheetCellType(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetCellTypeW';

function xlSheetReadStr(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): PWideChar cdecl;
external 'libxl' name 'xlSheetReadStrW';

function xlSheetWriteStr(handle: SheetHandle; row: Integer; col: Integer; value: PWideChar; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteStrW';

function xlSheetReadRichStr(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): RichStringHandle cdecl;
external 'libxl' name 'xlSheetReadRichStrW';

function xlSheetWriteRichStr(handle: SheetHandle; row: Integer; col: Integer; value: RichStringHandle; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteRichStrW';

function xlSheetReadNum(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): double cdecl;
external 'libxl' name 'xlSheetReadNumW';

function xlSheetWriteNum(handle: SheetHandle; row: Integer; col: Integer; value: double; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteNumW';

function xlSheetReadBool(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetReadBoolW';

function xlSheetWriteBool(handle: SheetHandle; row: Integer; col: Integer; value: Integer; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteBoolW';

function xlSheetReadBlank(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetReadBlankW';

function xlSheetWriteBlank(handle: SheetHandle; row: Integer; col: Integer; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteBlankW';

function xlSheetReadComment(handle: SheetHandle; row: Integer; col: Integer): PWideChar cdecl;
external 'libxl' name 'xlSheetReadCommentW';

procedure xlSheetWriteComment(handle: SheetHandle; row: Integer; col: Integer; const value: PWideChar; const author: PWideChar; width: Integer; height: Integer) cdecl;
external 'libxl' name 'xlSheetWriteCommentW';

procedure xlSheetRemoveComment(handle: SheetHandle; row: Integer; col: Integer) cdecl;
external 'libxl' name 'xlSheetRemoveCommentW';

function xlSheetIsDate(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetIsDateW';

function xlSheetIsRichStr(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetIsRichStrW';

function xlSheetReadError(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetReadErrorW';

procedure xlSheetWriteError(handle: SheetHandle; row: Integer; col: Integer; error: Integer; format: FormatHandle) cdecl;
external 'libxl' name 'xlSheetWriteErrorW';

function xlSheetColWidth(handle: SheetHandle; col: Integer): double cdecl;
external 'libxl' name 'xlSheetColWidthW';

function xlSheetRowHeight(handle: SheetHandle; row: Integer): double cdecl;
external 'libxl' name 'xlSheetRowHeightW';

function xlSheetSetCol(handle: SheetHandle; colFirst: Integer; colLast: Integer; width: double; format: FormatHandle; hidden: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetColW';

function xlSheetSetRow(handle: SheetHandle; row: Integer; height: double; format: FormatHandle; hidden: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetRowW';

function xlSheetRowHidden(handle: SheetHandle; row: Integer): Integer cdecl;
external 'libxl' name 'xlSheetRowHiddenW';

function xlSheetSetRowHidden(handle: SheetHandle; row: Integer; hidden: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetRowHiddenW';

function xlSheetColHidden(handle: SheetHandle; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetColHiddenW';

function xlSheetSetColHidden(handle: SheetHandle; col: Integer; hidden: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetColHiddenW';

function xlSheetGetMerge(handle: SheetHandle; row: Integer; col: Integer; rowFirst: PInteger; rowLast: PInteger; colFirst: PInteger; colLast: PInteger) : Integer cdecl;
external 'libxl' name 'xlSheetGetMergeW';

function xlSheetSetMerge(handle: SheetHandle; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetMergeW';

function xlSheetDelMerge(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetDelMergeW';

function xlSheetPictureSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetPictureSizeW';

function xlSheetGetPicture(handle: SheetHandle; index: Integer; rowTop: PInteger; colLeft: PInteger; rowBottom: PInteger; colRight: PInteger;
                                                 width: PInteger; height: PInteger; offset_x: PInteger; offset_y: PInteger): Integer cdecl;
external 'libxl' name 'xlSheetGetPictureW';

function xlSheetRemovePictureByIndex(handle: SheetHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlSheetRemovePictureByIndexW';

procedure xlSheetSetPicture(handle: SheetHandle; row: Integer; col: Integer; pictureId: Integer; scale: double; offset_x: Integer; offset_y: Integer; pos: Integer) cdecl;
external 'libxl' name 'xlSheetSetPictureW';

procedure xlSheetSetPicture2(handle: SheetHandle; row: Integer; col: Integer; pictureId: Integer; widht: Integer; height: Integer; offset_x: Integer; offset_y: Integer; pos: Integer) cdecl;
external 'libxl' name 'xlSheetSetPicture2W';

function xlSheetRemovePicture(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetRemovePictureW';

function xlSheetGetHorPageBreak(handle: SheetHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlSheetGetHorPageBreakW';

function xlSheetGetHorPageBreakSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetGetHorPageBreakSizeW';

function xlSheetGetVerPageBreak(handle: SheetHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlSheetGetVerPageBreakW';

function xlSheetGetVerPageBreakSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetGetVerPageBreakSizeW';

function xlSheetSetHorPageBreak(handle: SheetHandle; row: Integer; pageBreak: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetHorPageBreakW';

function xlSheetSetVerPageBreak(handle: SheetHandle; col: Integer; pageBreak: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetVerPageBreakW';

procedure xlSheetSplit(handle: SheetHandle; row: Integer; col: Integer) cdecl;
external 'libxl' name 'xlSheetSplitW';

function xlSheetGroupRows(handle: SheetHandle; rowFirst: Integer; rowLast: Integer; collapsed: Integer): Integer cdecl;
external 'libxl' name 'xlSheetGroupRowsW';

function xlSheetGroupCols(handle: SheetHandle; colFirst: Integer; colLast: Integer; collapsed: Integer): Integer cdecl;
external 'libxl' name 'xlSheetGroupColsW';

function xlSheetGroupSummaryBelow(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetGroupSummaryBelowW';

procedure xlSheetSetGroupSummaryBelow(handle: SheetHandle; below: Integer) cdecl;
external 'libxl' name 'xlSheetSetGroupSummaryBelowW';

function xlSheetGroupSummaryRight(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetGroupSummaryRightW';

procedure xlSheetSetGroupSummaryRight(handle: SheetHandle; right: Integer) cdecl;
external 'libxl' name 'xlSheetSetGroupSummaryRightW';

function xlSheetClear(handle: SheetHandle; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetClearW';

function xlSheetInsertCol(handle: SheetHandle; colFirst: Integer; colLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetInsertColW';

function xlSheetInsertRow(handle: SheetHandle; rowFirst: Integer; rowLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetInsertRowW';

function xlSheetRemoveCol(handle: SheetHandle; colFirst: Integer; colLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetRemoveColW';

function xlSheetRemoveRow(handle: SheetHandle; rowFirst: Integer; rowLast: Integer): Integer cdecl;
external 'libxl' name 'xlSheetRemoveRowW';

function xlSheetFirstRow(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetFirstRowW';

function xlSheetLastRow(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetLastRowW';

function xlSheetFirstCol(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetFirstColW';

function xlSheetLastCol(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetLastColW';

function xlSheetFirstFilledRow(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetFirstFilledRowW';

function xlSheetLastFilledRow(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetLastFilledRowW';

function xlSheetFirstFilledCol(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetFirstFilledColW';

function xlSheetLastFilledCol(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetLastFilledColW';

function xlSheetDisplayGridlines(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetDisplayGridlinesW';

procedure xlSheetSetDisplayGridlines(handle: SheetHandle; show: Integer) cdecl;
external 'libxl' name 'xlSheetSetDisplayGridlinesW';

function xlSheetPrintGridlines(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetPrintGridlinesW';

procedure xlSheetSetPrintGridlines(handle: SheetHandle; print: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintGridlinesW';

function xlSheetZoom(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetZoomW';

procedure xlSheetSetZoom(handle: SheetHandle; print: Integer) cdecl;
external 'libxl' name 'xlSheetSetZoomW';

function xlSheetPrintZoom(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetPrintZoomW';

procedure xlSheetSetPrintZoom(handle: SheetHandle; print: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintZoomW';

function xlSheetGetPrintFit(handle: SheetHandle; wPages: PInteger; hPages: PInteger): Integer cdecl;
external 'libxl' name 'xlSheetGetPrintFitW';

procedure xlSheetSetPrintFit(handle: SheetHandle; wPages: Integer; hPages: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintFitW';

function xlSheetLandscape(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetLandscapeW';

procedure xlSheetSetLandscape(handle: SheetHandle; landscape: Integer) cdecl;
external 'libxl' name 'xlSheetSetLandscapeW';

function xlSheetPaper(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetPaperW';

procedure xlSheetSetPaper(handle: SheetHandle; Paper: Integer) cdecl;
external 'libxl' name 'xlSheetSetPaperW';

function xlSheetHeader(handle: SheetHandle): PWideChar cdecl;
external 'libxl' name 'xlSheetHeaderW';

function xlSheetSetHeader(handle: SheetHandle; const header: PWideChar; margin: double): Integer cdecl;
external 'libxl' name 'xlSheetSetHeaderW';

function xlSheetHeaderMargin(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetHeaderMarginW';

function xlSheetFooter(handle: SheetHandle): PWideChar cdecl;
external 'libxl' name 'xlSheetFooterW';

function xlSheetSetFooter(handle: SheetHandle; const footer: PWideChar; margin: double): Integer cdecl;
external 'libxl' name 'xlSheetSetFooterW';

function xlSheetFooterMargin(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetFooterMarginW';

function xlSheetHCenter(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetHCenterW';

procedure xlSheetSetHCenter(handle: SheetHandle; hCenter: Integer) cdecl;
external 'libxl' name 'xlSheetSetHCenterW';

function xlSheetVCenter(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetVCenterW';

procedure xlSheetSetVCenter(handle: SheetHandle; vCenter: Integer) cdecl;
external 'libxl' name 'xlSheetSetVCenterW';

function xlSheetMarginLeft(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetMarginLeftW';

procedure xlSheetSetMarginLeft(handle: SheetHandle; margin: double) cdecl;
external 'libxl' name 'xlSheetSetMarginLeftW';

function xlSheetMarginRight(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetMarginRightW';

procedure xlSheetSetMarginRight(handle: SheetHandle; margin: double) cdecl;
external 'libxl' name 'xlSheetSetMarginRightW';

function xlSheetMarginTop(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetMarginTopW';

procedure xlSheetSetMarginTop(handle: SheetHandle; margin: double) cdecl;
external 'libxl' name 'xlSheetSetMarginTopW';

function xlSheetMarginBottom(handle: SheetHandle): double cdecl;
external 'libxl' name 'xlSheetMarginBottomW';

procedure xlSheetSetMarginBottom(handle: SheetHandle; margin: double) cdecl;
external 'libxl' name 'xlSheetSetMarginBottomW';

function xlSheetPrintRowCol(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetPrintRowColW';

procedure xlSheetSetPrintRowCol(handle: SheetHandle; print: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintRowColW';

function xlSheetProtect(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetProtectW';

procedure xlSheetSetProtect(handle: SheetHandle; protect: Integer) cdecl;
external 'libxl' name 'xlSheetSetProtectW';

procedure xlSheetSetProtectEx(handle: SheetHandle; protect: Integer; const name: PWideChar; enhancedProtection: Integer) cdecl;
external 'libxl' name 'xlSheetSetProtectExW';

function xlSheetRightToLeft(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetRightToLeftW';

procedure xlSheetSetRightToLeft(handle: SheetHandle; rightToLeft: Integer) cdecl;
external 'libxl' name 'xlSheetSetRightToLeftW';

procedure xlSheetSetPrintRepeatRows(handle: SheetHandle; rowFirst: Integer; rowLast: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintRepeatRowsW';

procedure xlSheetSetPrintRepeatCols(handle: SheetHandle; colFirst: Integer; colLast: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintRepeatColsW';

procedure xlSheetSetPrintArea(handle: SheetHandle; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer) cdecl;
external 'libxl' name 'xlSheetSetPrintAreaW';

procedure xlSheetClearPrintRepeats(handle: SheetHandle) cdecl;
external 'libxl' name 'xlSheetClearPrintRepeatsW';

procedure xlSheetClearPrintArea(handle: SheetHandle) cdecl;
external 'libxl' name 'xlSheetClearPrintAreaW';

function xlSheetGetNamedRange(handle: SheetHandle; const name: PWideChar; rowFirst: PInteger; rowLast: PInteger; colFirst: PInteger; colLast: PInteger; scopeId: Integer; hidden: PBoolean): Integer cdecl;
external 'libxl' name 'xlSheetGetNamedRangeW';

function xlSheetSetNamedRange(handle: SheetHandle; const name: PWideChar; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; scopeId: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetNamedRangeW';

function xlSheetDelNamedRange(handle: SheetHandle; const name: PWideChar; scopeId: Integer): Integer cdecl;
external 'libxl' name 'xlSheetDelNamedRangeW';

function xlSheetNamedRangeSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetNamedRangeSizeW';

function xlSheetNamedRange(handle: SheetHandle; index: Integer; rowFirst: PInteger; rowLast: PInteger; colFirst: PInteger; colLast: PInteger; scopeId: PInteger; hidden: PBoolean): PWideChar cdecl;
external 'libxl' name 'xlSheetNamedRangeW';

function xlSheetTableSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetTableSizeW';

function xlSheetTable(handle: SheetHandle; index: Integer; rowFirst: PInteger; rowLast: PInteger; colFirst: PInteger; colLast: PInteger; headerRowCount: PInteger; totalsRowCount: PInteger): PWideChar cdecl;
external 'libxl' name 'xlSheetTableW';

function xlSheetHyperlinkSize(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetHyperlinkSizeW';

function xlSheetHyperlink(handle: SheetHandle; index: Integer; rowFirst: PInteger; rowLast: PInteger; colFirst: PInteger; colLast: PInteger): PWideChar cdecl;
external 'libxl' name 'xlSheetHyperlinkW';

function xlSheetDelHyperlink(handle: SheetHandle; index: Integer): Integer cdecl;
external 'libxl' name 'xlSheetDelHyperlinkW';

procedure xlSheetAddHyperlink(handle: SheetHandle; const hyperlink: PWideChar; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer) cdecl;
external 'libxl' name 'xlSheetAddHyperlinkW';

function xlSheetAutoFilter(handle: SheetHandle): AutoFilterHandle cdecl;
external 'libxl' name 'xlSheetAutoFilterW';

procedure xlSheetApplyFilter(handle: SheetHandle) cdecl;
external 'libxl' name 'xlSheetApplyFilterW';

procedure xlSheetRemoveFilter(handle: SheetHandle) cdecl;
external 'libxl' name 'xlSheetRemoveFilterW';

function xlSheetName(handle: SheetHandle): PWideChar cdecl;
external 'libxl' name 'xlSheetNameW';

procedure xlSheetSetName(handle: SheetHandle; const name: PWideChar) cdecl;
external 'libxl' name 'xlSheetSetNameW';

function xlSheetCellFormat(handle: SheetHandle; row: Integer; col: Integer): FormatHandle cdecl;
external 'libxl' name 'xlSheetCellFormatW';

procedure xlSheetSetCellFormat(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle) cdecl;
external 'libxl' name 'xlSheetSetCellFormatW';

function xlSheetIsFormula(handle: SheetHandle; row: Integer; col: Integer): Integer cdecl;
external 'libxl' name 'xlSheetIsFormulaW';

function xlSheetReadFormula(handle: SheetHandle; row: Integer; col: Integer; format: PFormatHandle): PWideChar cdecl;
external 'libxl' name 'xlSheetReadFormulaW';

function xlSheetWriteFormula(handle: SheetHandle; row: Integer; col: Integer; value: PWideChar; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteFormulaW';

function xlSheetWriteFormulaNum(handle: SheetHandle; row: Integer; col: Integer; expr: PWideChar; value: double; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteFormulaNumW';

function xlSheetWriteFormulaStr(handle: SheetHandle; row: Integer; col: Integer; expr: PWideChar; value: PWideChar; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteFormulaStrW';

function xlSheetWriteFormulaBool(handle: SheetHandle; row: Integer; col: Integer; expr: PWideChar; value: Integer; format: FormatHandle): Integer cdecl;
external 'libxl' name 'xlSheetWriteFormulaBoolW';

function xlSheetHidden(handle: SheetHandle): Integer cdecl;
external 'libxl' name 'xlSheetHiddenW';

function xlSheetSetHidden(handle: SheetHandle; hidden: Integer): Integer cdecl;
external 'libxl' name 'xlSheetSetHiddenW';

procedure xlSheetGetTopLeftView(handle: SheetHandle; row: PInteger; col: PInteger) cdecl;
external 'libxl' name 'xlSheetGetTopLeftViewW';

procedure xlSheetSetTopLeftView(handle: SheetHandle; row: Integer; col: Integer) cdecl;
external 'libxl' name 'xlSheetSetTopLeftViewW';

procedure xlSheetSetAutoFitArea(handle: SheetHandle; rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer) cdecl;
external 'libxl' name 'xlSheetSetAutoFitAreaW';

procedure xlSheetAddrToRowCol(handle: SheetHandle; const name: PWideChar; row: PInteger; col: PInteger; rowRelative: PInteger; colRelative: PInteger) cdecl;
external 'libxl' name 'xlSheetAddrToRowColW';

function xlSheetRowColToAddr(handle: SheetHandle; row: Integer; col: Integer; rowRelative: Integer; colRelative: Integer): PWideChar cdecl;
external 'libxl' name 'xlSheetRowColToAddrW';

procedure xlSheetSetTabColor(handle: SheetHandle; color: Integer) cdecl;
external 'libxl' name 'xlSheetSetTabColorW';

procedure xlSheetSetTabRgbColor(handle: SheetHandle; red: Integer; green: Integer; blue: Integer) cdecl;
external 'libxl' name 'xlSheetSetTabRgbColorW';

function xlSheetAddIgnoredError(handle: SheetHandle; rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer; iError: Integer): Integer cdecl;
external 'libxl' name 'xlSheetAddIgnoredErrorW';

procedure xlSheetAddDataValidation(handle: SheetHandle; vtype: Integer; op: Integer; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: PWideChar; value2: PWideChar) cdecl;
external 'libxl' name 'xlSheetAddDataValidationW';

procedure xlSheetAddDataValidationEx(handle: SheetHandle; vtype: Integer; op: Integer; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: PWideChar; value2: PWideChar;
                                     allowBlank: Integer; hideDropDown: Integer; showInputMessage: Integer; showErrorMessage: Integer; promptTitle: PWideChar; prompt: PWideChar;
                                     errorTitle: PWideChar; error: PWideChar; errorStyle: Integer) cdecl;
external 'libxl' name 'xlSheetAddDataValidationExW';

procedure xlSheetAddDataValidationDouble(handle: SheetHandle; vtype: Integer; op: Integer; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double) cdecl;
external 'libxl' name 'xlSheetAddDataValidationDoubleW';

procedure xlSheetAddDataValidationDoubleEx(handle: SheetHandle; vtype: Integer; op: Integer; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double;
                                     allowBlank: Integer; hideDropDown: Integer; showInputMessage: Integer; showErrorMessage: Integer; promptTitle: PWideChar; prompt: PWideChar;
                                     errorTitle: PWideChar; error: PWideChar; errorStyle: Integer) cdecl;
external 'libxl' name 'xlSheetAddDataValidationDoubleExW';

procedure xlSheetRemoveDataValidations(handle: SheetHandle) cdecl;
external 'libxl' name 'xlSheetRemoveDataValidationsW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlFormatFont(handle: FormatHandle): FontHandle cdecl;
external 'libxl' name 'xlFormatFontW';

function xlFormatSetFont(handle: FormatHandle; FontHandle: FontHandle): FontHandle cdecl;
external 'libxl' name 'xlFormatSetFontW';

function xlFormatNumFormat(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatNumFormatW';

procedure xlFormatSetNumFormat(handle: FormatHandle; NumFormat: Integer) cdecl;
external 'libxl' name 'xlFormatSetNumFormatW';

function xlFormatAlignH(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatAlignHW';

procedure xlFormatSetAlignH(handle: FormatHandle; align: Integer) cdecl;
external 'libxl' name 'xlFormatSetAlignHW';

function xlFormatAlignV(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatAlignVW';

procedure xlFormatSetAlignV(handle: FormatHandle; align: Integer) cdecl;
external 'libxl' name 'xlFormatSetAlignVW';

function xlFormatWrap(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatWrapW';

procedure xlFormatSetWrap(handle: FormatHandle; wrap: Integer) cdecl;
external 'libxl' name 'xlFormatSetWrapW';

function xlFormatRotation(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatRotationW';

function xlFormatSetRotation(handle: FormatHandle; rotation: Integer): Integer cdecl;
external 'libxl' name 'xlFormatSetRotationW';

function xlFormatIndent(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatIndentW';

procedure xlFormatSetIndent(handle: FormatHandle; indent: Integer) cdecl;
external 'libxl' name 'xlFormatSetIndentW';

function xlFormatShrinkToFit(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatShrinkToFitW';

procedure xlFormatSetShrinkToFit(handle: FormatHandle; shrinkToFit: Integer) cdecl;
external 'libxl' name 'xlFormatSetShrinkToFitW';

procedure xlFormatSetBorder(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderW';

procedure xlFormatSetBorderColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderColorW';

function xlFormatBorderLeft(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderLeftW';

procedure xlFormatSetBorderLeft(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderLeftW';

function xlFormatBorderRight(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderRightW';

procedure xlFormatSetBorderRight(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderRightW';

function xlFormatBorderTop(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderTopW';

procedure xlFormatSetBorderTop(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderTopW';

function xlFormatBorderBottom(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderBottomW';

procedure xlFormatSetBorderBottom(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderBottomW';

function xlFormatBorderLeftColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderLeftColorW';

procedure xlFormatSetBorderLeftColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderLeftColorW';

function xlFormatBorderRightColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderRightColorW';

procedure xlFormatSetBorderRightColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderRightColorW';

function xlFormatBorderTopColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderTopColorW';

procedure xlFormatSetBorderTopColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderTopColorW';

function xlFormatBorderBottomColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderBottomColorW';

procedure xlFormatSetBorderBottomColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderBottomColorW';

function xlFormatBorderDiagonal(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderDiagonalW';

procedure xlFormatSetBorderDiagonal(handle: FormatHandle; border: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderDiagonalW';

function xlFormatBorderDiagonalStyle(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderDiagonalStyleW';

procedure xlFormatSetBorderDiagonalStyle(handle: FormatHandle; style: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderDiagonalStyleW';

function xlFormatBorderDiagonalColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatBorderDiagonalColorW';

procedure xlFormatSetBorderDiagonalColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetBorderDiagonalColorW';

function xlFormatFillPattern(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatFillPatternW';

procedure xlFormatSetFillPattern(handle: FormatHandle; pattern: Integer) cdecl;
external 'libxl' name 'xlFormatSetFillPatternW';

function xlFormatPatternForegroundColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatPatternForegroundColorW';

procedure xlFormatSetPatternForegroundColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetPatternForegroundColorW';

function xlFormatPatternBackgroundColor(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatPatternBackgroundColorW';

procedure xlFormatSetPatternBackgroundColor(handle: FormatHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFormatSetPatternBackgroundColorW';

function xlFormatLocked(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatLockedW';

procedure xlFormatSetLocked(handle: FormatHandle; locked: Integer) cdecl;
external 'libxl' name 'xlFormatSetLockedW';

function xlFormatHidden(handle: FormatHandle): Integer cdecl;
external 'libxl' name 'xlFormatHiddenW';

procedure xlFormatSetHidden(handle: FormatHandle; hidden: Integer) cdecl;
external 'libxl' name 'xlFormatSetHiddenW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlFontSize(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontSizeW';

procedure xlFontSetSize(handle: FontHandle; size: Integer) cdecl;
external 'libxl' name 'xlFontSetSizeW';

function xlFontItalic(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontItalicW';

procedure xlFontSetItalic(handle: FontHandle; italic: Integer) cdecl;
external 'libxl' name 'xlFontSetItalicW';

function xlFontStrikeOut(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontStrikeOutW';

procedure xlFontSetStrikeOut(handle: FontHandle; strikeOut: Integer) cdecl;
external 'libxl' name 'xlFontSetStrikeOutW';

function xlFontColor(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontColorW';

procedure xlFontSetColor(handle: FontHandle; Color: Integer) cdecl;
external 'libxl' name 'xlFontSetColorW';

function xlFontBold(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontBoldW';

procedure xlFontSetBold(handle: FontHandle; bold: Integer) cdecl;
external 'libxl' name 'xlFontSetBoldW';

function xlFontScript(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontScriptW';

procedure xlFontSetScript(handle: FontHandle; Script: Integer) cdecl;
external 'libxl' name 'xlFontSetScriptW';

function xlFontUnderline(handle: FontHandle): Integer cdecl;
external 'libxl' name 'xlFontUnderlineW';

procedure xlFontSetUnderline(handle: FontHandle; Underline: Integer) cdecl;
external 'libxl' name 'xlFontSetUnderlineW';

function xlFontName(handle: FontHandle): PWideChar cdecl;
external 'libxl' name 'xlFontNameW';

function xlFontSetName(handle: FontHandle; const name: PWideChar): Integer cdecl;
external 'libxl' name 'xlFontSetNameW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlAutoFilterGetRef(handle: AutoFilterHandle; rowFirst, rowLast, colFirst, colLast: PInteger): Integer cdecl;
external 'libxl' name 'xlAutoFilterGetRefW';

procedure xlAutoFilterSetRef(handle: AutoFilterHandle; rowFirst, rowLast, colFirst, colLast: Integer) cdecl;
external 'libxl' name 'xlAutoFilterSetRefW';

function xlAutoFilterColumn(handle: AutoFilterHandle; colId: Integer): FilterColumnHandle cdecl;
external 'libxl' name 'xlAutoFilterColumnW';

function xlAutoFilterColumnSize(handle: AutoFilterHandle): Integer cdecl;
external 'libxl' name 'xlAutoFilterColumnSizeW';

function xlAutoFilterColumnByIndex(handle: AutoFilterHandle; index: Integer): FilterColumnHandle cdecl;
external 'libxl' name 'xlAutoFilterColumnByIndexW';

function xlAutoFilterGetSortRange(handle: AutoFilterHandle; rowFirst, rowLast, colFirst, colLast: PInteger): Integer cdecl;
external 'libxl' name 'xlAutoFilterGetSortRangeW';

function xlAutoFilterGetSort(handle: AutoFilterHandle; columnIndex, descending: PInteger): Integer cdecl;
external 'libxl' name 'xlAutoFilterGetSortW';

function xlAutoFilterSetSort(handle: AutoFilterHandle; columnIndex, descending: Integer): Integer cdecl;
external 'libxl' name 'xlAutoFilterSetSortW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlFilterColumnIndex(handle: FilterColumnHandle): Integer cdecl;
external 'libxl' name 'xlFilterColumnIndexW';

function xlFilterColumnFilterType(handle: FilterColumnHandle): Integer cdecl;
external 'libxl' name 'xlFilterColumnFilterTypeW';

function xlFilterColumnFilterSize(handle: FilterColumnHandle): Integer cdecl;
external 'libxl' name 'xlFilterColumnFilterSizeW';

function xlFilterColumnFilter(handle: FilterColumnHandle; index: Integer): PWideChar cdecl;
external 'libxl' name 'xlFilterColumnFilterW';

procedure xlFilterColumnAddFilter(handle: FilterColumnHandle; value: PWideChar) cdecl;
external 'libxl' name 'xlFilterColumnAddFilterW';

function xlFilterColumnGetTop10(handle: FilterColumnHandle; value: PDouble; top, percent: PInteger): Integer cdecl;
external 'libxl' name 'xlFilterColumnGetTop10W';

procedure xlFilterColumnSetTop10(handle: FilterColumnHandle; value: double; top, percent: Integer) cdecl;
external 'libxl' name 'xlFilterColumnSetTop10W';

function xlFilterColumnGetCustomFilter(handle: FilterColumnHandle; op1: PInteger; v1: PAnsiString; op2: PInteger; v2: PAnsiString; andOp: PInteger): Integer cdecl;
external 'libxl' name 'xlFilterColumnGetCustomFilterW';

procedure xlFilterColumnSetCustomFilter(handle: FilterColumnHandle; op: Integer; v: WideString) cdecl;
external 'libxl' name 'xlFilterColumnSetCustomFilterW';

procedure xlFilterColumnSetCustomFilterEx(handle: FilterColumnHandle; op1: Integer; v1: WideString; op2: Integer; v2: WideString; andOp: Integer) cdecl;
external 'libxl' name 'xlFilterColumnSetCustomFilterExW';

procedure xlFilterColumnClear(handle: FilterColumnHandle) cdecl;
external 'libxl' name 'xlFilterColumnClearW';

////////////////////////////////////////////////////////////////////////////////////////////////////

function xlRichStringAddFont(handle: RichStringHandle; initFont: FontHandle): FontHandle cdecl;
external 'libxl' name 'xlRichStringAddFontW';

procedure xlRichStringAddText(handle: RichStringHandle; text: PWideChar; font: FontHandle) cdecl;
external 'libxl' name 'xlRichStringAddTextW';

function xlRichStringGetText(handle: RichStringHandle; index: Integer; font: PFontHandle): PWideChar cdecl;
external 'libxl' name 'xlRichStringGetTextW';

function xlRichStringTextSize(handle: RichStringHandle): Integer cdecl;
external 'libxl' name 'xlRichStringTextSizeW';

////////////////////////////////////////////////////////////////////////////////////////////////////


{$M+}

type
  TXLFilterColumn = class(TObject)
  private
    book: BookHandle;

  protected
  public
     handle: FilterColumnHandle;

     function index(): Integer;

     function filterType(): Filter;

     function filterSize(): Integer;
     function filter(index: Integer): PWideChar;
     procedure addFilter(const value: PWideChar);

     function getTop10(var value: double; var top, percent: boolean): boolean;
     procedure setTop10(value: double); overload;
     procedure setTop10(value: double; top: boolean); overload;
     procedure setTop10(value: double; top, percent: boolean); overload;

     function getCustomFilter(var op1: Operator; var v1: PWideChar; var op2: Operator; var v2: PWideChar; var andOp: boolean): boolean;
     procedure setCustomFilter(op: Operator; const v: PWideChar); overload;
     procedure setCustomFilter(op1: Operator; const v1: PWideChar; op2: Operator; const v2: PWideChar; andOp: boolean); overload;

     procedure clear();

  published
    constructor Create(handle: FilterColumnHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLAutoFilter = class(TObject)
  private
    book: BookHandle;
    filterColumns: array of TXLFilterColumn;

  protected
  public
     handle: AutoFilterHandle;

     function getRef(var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
     procedure setRef(rowFirst, rowLast, colFirst, colLast: Integer);

     function column(colId: Integer): TXLFilterColumn;

     function columnSize(): Integer;
     function columnByIndex(index: Integer): TXLFilterColumn;

     function getSortRange(var rowFirst, rowLast, colFirst, colLast: Integer): boolean;

     function getSort(var columnIndex: Integer; var descending: boolean): boolean;
     function setSort(columnIndex: Integer; descending: boolean): boolean;

  published
    constructor Create(handle: AutoFilterHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLFont = class(TObject)
  private
    book: BookHandle;

    procedure Setsize(const value: Integer);
    function Getsize: Integer;

    procedure Setbold(const value: boolean);
    function Getbold: boolean;

    procedure Setitalic(const value: boolean);
    function Getitalic: boolean;

    procedure SetUnderline(const value: Underline);
    function GetUnderline: Underline;

    procedure SetstrikeOut(const value: boolean);
    function GetstrikeOut: boolean;

    procedure Setscript(const value: Script);
    function Getscript: Script;

    procedure Setcolor(const value: Color);
    function Getcolor: Color;

    procedure SetName(const value: WideString);
    function GetName: WideString;

  protected
  public
    handle: FontHandle;

    property size: Integer read Getsize write Setsize;
    property bold: boolean read Getbold write Setbold;
    property italic: boolean read Getitalic write Setitalic;
    property underline: Underline read GetUnderline write SetUnderline;
    property strikeOut: boolean read GetstrikeOut write SetstrikeOut;
    property script: Script read Getscript write Setscript;
    property color: Color read Getcolor write Setcolor;
    property name: WideString read GetName write SetName;

  published
    constructor Create(handle: FontHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLRichString = class(TObject)
  private
    book: BookHandle;
    fonts: array of TXLFont;

  protected
  public
     handle: RichStringHandle;

     function addFont(): TXLFont; overload;
     function addFont(initFont: TXLFont): TXLFont; overload;
     procedure addText(const text: PWideChar); overload;
     procedure addText(const text: PWideChar; font: TXLFont); overload;
     function getText(index: Integer): WideString; overload;
     function getText(index: Integer; var font: TXLFont): WideString; overload;
     function textSize(): Integer;

  published
    constructor Create(handle: RichStringHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLFormat = class(TObject)
  private
    book: BookHandle;
    ffont: TXLFont;

    procedure Setfont(const value: TXLFont);
    function Getfont: TXLFont;

    procedure SetalignH(const value: AlignH);
    function GetalignH: AlignH;

    procedure SetAlignV(const value: AlignV);
    function GetAlignV: AlignV;
    procedure Setwrap(const value: boolean);
    function Getwrap: boolean;
    procedure Setrotation(const value: Integer);
    function Getrotation: Integer;
    procedure Setindent(const value: Integer);
    function Getindent: Integer;
    procedure SetshrinkToFit(const value: boolean);
    function GetshrinkToFit: boolean;

    procedure SetborderLeft(const value: BorderStyle);
    function GetborderLeft: BorderStyle;
    procedure SetborderRight(const value: BorderStyle);
    function GetborderRight: BorderStyle;
    procedure SetborderTop(const value: BorderStyle);
    function GetborderTop: BorderStyle;
    procedure SetborderBottom(const value: BorderStyle);
    function GetborderBottom: BorderStyle;

    procedure SetborderLeftColor(const value: Color);
    function GetborderLeftColor: Color;
    procedure SetborderRightColor(const value: Color);
    function GetborderRightColor: Color;
    procedure SetborderTopColor(const value: Color);
    function GetborderTopColor: Color;
    procedure SetborderBottomColor(const value: Color);
    function GetborderBottomColor: Color;

    procedure SetborderDiagonal(const value: BorderDiagonal);
    function GetborderDiagonal: BorderDiagonal;
    procedure SetBorderDiagonalStyle(const value: BorderStyle);
    function GetBorderDiagonalStyle: BorderStyle;
    procedure SetborderDiagonalColor(const value: Color);
    function GetborderDiagonalColor: Color;

    procedure SetfillPattern(const value: FillPattern);
    function GetfillPattern: FillPattern;

    procedure SetpatternForegroundColor(const value: Color);
    function GetpatternForegroundColor: Color;
    procedure SetpatternBackgroundColor(const value: Color);
    function GetpatternBackgroundColor: Color;

    procedure Setlocked(const value: boolean);
    function Getlocked: boolean;

    procedure Sethidden(const value: boolean);
    function Gethidden: boolean;

  protected
  public
    handle: FormatHandle;

    function numFormat: Integer;
    procedure setNumFormat(NumFormat: Integer); overload;
    procedure setNumFormat(value: NumFormat); overload;
    procedure setBorder(style: BorderStyle);
    procedure setBorderColor(value: Color);

    property hidden: boolean read Gethidden write Sethidden;
    property locked: boolean read Getlocked write Setlocked;
    property fillPattern: FillPattern read GetfillPattern write SetfillPattern;
    property patternForegroundColor: Color read GetpatternForegroundColor write SetpatternForegroundColor;
    property patternBackgroundColor: Color read GetpatternBackgroundColor write SetpatternBackgroundColor;
    property borderDiagonal: BorderDiagonal read GetborderDiagonal write SetborderDiagonal;
    property borderDiagonalStyle: BorderStyle read GetBorderDiagonalStyle write SetBorderDiagonalStyle;
    property borderDiagonalColor: Color read GetborderDiagonalColor write SetborderDiagonalColor;

    property font: TXLFont read Getfont write Setfont;
    property alignH: AlignH read GetalignH write SetalignH;
    property alignV: AlignV read GetAlignV write SetAlignV;
    property wrap: boolean read Getwrap write Setwrap;
    property rotation: Integer read Getrotation write Setrotation;
    property indent: Integer read Getindent write Setindent;
    property shrinkToFit: boolean read GetshrinkToFit write SetshrinkToFit;

    property borderLeft: BorderStyle read GetborderLeft write SetborderLeft;
    property borderRight: BorderStyle read GetborderRight write SetborderRight;
    property borderTop: BorderStyle read GetborderTop write SetborderTop;
    property borderBottom: BorderStyle read GetborderBottom write SetborderBottom;

    property borderLeftColor: Color read GetborderLeftColor write SetborderLeftColor;
    property borderRightColor: Color read GetborderRightColor write SetborderRightColor;
    property borderTopColor: Color read GetborderTopColor write SetborderTopColor;
    property borderBottomColor: Color read GetborderBottomColor write SetborderBottomColor;

  published
    constructor Create(handle: FormatHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLSheet = class(TObject)
  private
    handle: SheetHandle;
    book: BookHandle;
    richStrings: array of TXLRichString;
    formats: array of TXLFormat;
    xlAutoFilter: TXLAutoFilter;

    procedure SetDisplayGridlines(const value: boolean);
    function GetDisplayGridlines: boolean;

    procedure SetprintGridlines(const value: boolean);
    function GetprintGridlines: boolean;

    procedure SetZoom(const value: Integer);
    function GetZoom: Integer;

    procedure SetprintZoom(const value: Integer);
    function GetprintZoom: Integer;

    procedure Setlandscape(const value: boolean);
    function Getlandscape: boolean;

    procedure SetPaper(const value: Paper);
    function GetPaper: Paper;

    procedure Setheader(const value: WideString);
    function Getheader: WideString;

    procedure SetheaderMargin(const value: double);
    function GetheaderMargin: double;

    function Getfooter: WideString;
    procedure Setfooter(const value: WideString);
    function GetfooterMargin: double;
    procedure SetfooterMargin(const value: double);

    function GethCenter: boolean;
    procedure SethCenter(const value: boolean);

    function GetvCenter: boolean;
    procedure SetvCenter(const value: boolean);

    function GetmarginLeft: double;
    procedure SetmarginLeft(const value: double);
    function GetmarginRight: double;
    procedure SetmarginRight(const value: double);
    function GetmarginBottom: double;
    function GetmarginTop: double;
    procedure SetmarginBottom(const value: double);
    procedure SetmarginTop(const value: double);

    function GetprintRowCol: boolean;
    procedure SetprintRowCol(const value: boolean);

    function GetprotectImpl: boolean;
    procedure SetprotectImpl(const value: boolean);

    function Getrighttoleft: boolean;
    procedure Setrighttoleft(const value: boolean);

    function GetHidden: SheetState;
    procedure SetHidden(const value: SheetState);

    function GetName: WideString;
    procedure SetName(const value: WideString);

    function GetGroupSummaryBelow: boolean;
    procedure SetGroupSummaryBelow(const value: boolean);

    function GetGroupSummaryRight: boolean;
    procedure SetGroupSummaryRight(const value: boolean);

  protected
  public

    property displayGridlines: boolean read GetDisplayGridlines write SetDisplayGridlines;
    property printGridlines: boolean read GetprintGridlines write SetprintGridlines;
    property zoom: Integer read GetZoom write SetZoom;
    property printZoom: Integer read GetprintZoom write SetprintZoom;
    property landscape: boolean read Getlandscape write Setlandscape;
    property paper: Paper read GetPaper write SetPaper;
    property header: WideString read Getheader write Setheader;
    property headerMargin: double read GetheaderMargin write SetheaderMargin;
    property footer: WideString read Getfooter write Setfooter;
    property footerMargin: double read GetfooterMargin write SetfooterMargin;
    property hCenter: boolean read GethCenter write SethCenter;
    property vCenter: boolean read GetvCenter write SetvCenter;

    property marginLeft: double read GetmarginLeft write SetmarginLeft;
    property marginRight: double read GetmarginRight write SetmarginRight;
    property marginTop: double read GetmarginTop write SetmarginTop;
    property marginBottom: double read GetmarginBottom write SetmarginBottom;

    property printRowCol: boolean read GetprintRowCol write SetprintRowCol;
    property name: WideString read GetName write SetName;
    property groupSummaryBelow: boolean read GetGroupSummaryBelow write SetGroupSummaryBelow;
    property groupSummaryRight: boolean read GetGroupSummaryRight write SetGroupSummaryRight;

    property protect: boolean read GetprotectImpl write SetprotectImpl;
    property hidden: SheetState read GetHidden write SetHidden;
    property rightToLeft: boolean read Getrighttoleft write Setrighttoleft;

    function errorMessage(): WideString;

    function getCellType(row: Integer; col: Integer): CellType;
    function isFormula(row: Integer; col: Integer): boolean;

    function cellFormat(row: Integer; col: Integer): TXLFormat;
    procedure setCellFormat(row: Integer; col: Integer; format: TXLFormat);
    
    function readStr(row: Integer; col: Integer): WideString; overload;
    function readStr(row: Integer; col: Integer; var format: TXLFormat): WideString; overload;
    function writeStr(row: Integer; col: Integer; value: PWideChar): boolean; overload;
    function writeStr(row: Integer; col: Integer; value: PWideChar; format: TXLFormat): boolean; overload;

    function readRichStr(row: Integer; col: Integer): TXLRichString; overload;
    function readRichStr(row: Integer; col: Integer; var format: TXLFormat): TXLRichString; overload;
    function writeRichStr(row: Integer; col: Integer; value: TXLRichString): boolean; overload;
    function writeRichStr(row: Integer; col: Integer; value: TXLRichString; format: TXLFormat): boolean; overload;

    function readNum(row: Integer; col: Integer): double; overload;
    function readNum(row: Integer; col: Integer; var format: TXLFormat): double; overload;
    function writeNum(row: Integer; col: Integer; value: double): boolean; overload;
    function writeNum(row: Integer; col: Integer; value: double; format: TXLFormat): boolean; overload;

    function readBool(row: Integer; col: Integer): boolean; overload;
    function readBool(row: Integer; col: Integer; var format: TXLFormat): boolean; overload;
    function writeBool(row: Integer; col: Integer; value: boolean): boolean; overload;
    function writeBool(row: Integer; col: Integer; value: boolean; format: TXLFormat): boolean; overload;

    function readBlank(row: Integer; col: Integer; var format: TXLFormat): boolean;
    function writeBlank(row: Integer; col: Integer; format: TXLFormat): boolean;

    function readFormula(row: Integer; col: Integer): WideString; overload;
    function readFormula(row: Integer; col: Integer; var format: TXLFormat): WideString; overload;
    function writeFormula(row: Integer; col: Integer; value: PWideChar): boolean; overload;
    function writeFormula(row: Integer; col: Integer; value: PWideChar; format: TXLFormat): boolean; overload;

    function writeFormulaNum(row: Integer; col: Integer; expr: PWideChar; value: double): boolean; overload;
    function writeFormulaNum(row: Integer; col: Integer; expr: PWideChar; value: double; format: TXLFormat): boolean; overload;

    function writeFormulaStr(row: Integer; col: Integer; expr: PWideChar; value: PWideChar): boolean; overload;
    function writeFormulaStr(row: Integer; col: Integer; expr: PWideChar; value: PWideChar; format: TXLFormat): boolean; overload;

    function writeFormulaBool(row: Integer; col: Integer; expr: PWideChar; value: Integer): boolean; overload;
    function writeFormulaBool(row: Integer; col: Integer; expr: PWideChar; value: Integer; format: TXLFormat): boolean; overload;

    function readComment(row: Integer; col: Integer): WideString;
    procedure writeComment(row: Integer; col: Integer; value: PWideChar; author: PWideChar; width, height: Integer);
    procedure removeComment(row: Integer; col: Integer);

    function isDate(row, col: Integer): boolean;
    function isRichStr(row, col: Integer): boolean;

    function readError(row, col: Integer): ErrorType;
    procedure writeError(row, col: Integer; error: ErrorType); overload;
    procedure writeError(row, col: Integer; error: ErrorType; format: TXLFormat); overload;

    function colWidth(col: Integer): double;
    function rowHeight(col: Integer): double;

    function setCol(col: Integer; width: double): boolean; overload;
    function setCol(colFirst, colLast: Integer; width: double): boolean; overload;
    function setCol(col: Integer; width: double; format: TXLFormat): boolean; overload;
    function setCol(colFirst, colLast: Integer; width: double; format: TXLFormat): boolean; overload;
    function setCol(col: Integer; width: double; format: TXLFormat; hidden: boolean): boolean; overload;
    function setCol(colFirst, colLast: Integer; width: double; format: TXLFormat; hidden: boolean): boolean; overload;

    function setRow(row: Integer; height: double): boolean; overload;
    function setRow(row: Integer; height: double; format: TXLFormat): boolean; overload;
    function setRow(row: Integer; height: double; format: TXLFormat; hidden: boolean): boolean; overload;

    function rowHidden(row: Integer): boolean;
    function setRowHidden(row: Integer; hidden: boolean): boolean;

    function colHidden(col: Integer): boolean;
    function setColHidden(col: Integer; hidden: boolean): boolean;

    function getMerge(row, col: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
    function setMerge(rowFirst, rowLast, colFirst, colLast: Integer): boolean;
    function delMerge(row, col: Integer): boolean;

    function pictureSize: Integer;
    function getPicture(index: Integer): Integer; overload;
    function getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight: Integer): Integer; overload;
    function getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight, width, height: Integer): Integer; overload;
    function getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight, width, height, offset_x, offset_y: Integer): Integer; overload;
    function removePictureByIndex(index: Integer): boolean;

    procedure setPicture(row, col, pictureId: Integer); overload;
    procedure setPicture(row, col, pictureId: Integer; scale: double); overload;
    procedure setPicture(row, col, pictureId, width, height: Integer); overload;
    function removePicture(row, col: Integer): boolean;

    function getHorPageBreak(index: Integer): Integer;
    function getHorPageBreakSize: Integer;

    function getVerPageBreak(index: Integer): Integer;
    function getVerPageBreakSize: Integer;

    function setHorPageBreak(row: Integer): boolean;
    function delHorPageBreak(row: Integer): boolean;
    function setVerPageBreak(col: Integer): boolean;
    function delVerPageBreak(col: Integer): boolean;

    procedure split(row, col: Integer);

    function groupRows(rowFirst, rowLast: Integer; collapsed: boolean): boolean;
    function groupCols(colFirst, colLast: Integer; collapsed: boolean): boolean;

    function clear(rowFirst, rowLast: Integer; colFirst, colLast: Integer): boolean;

    function insertCol(colFirst, colLast: Integer): boolean;
    function insertRow(rowFirst, rowLast: Integer): boolean;
    function removeCol(colFirst, colLast: Integer): boolean;
    function removeRow(rowFirst, rowLast: Integer): boolean;

    function copyCell(rowSrc, colSrc, rowDst, colDst: Integer): boolean;

    function firstRow: Integer;
    function lastRow: Integer;
    function firstCol: Integer;
    function lastCol: Integer;

    function firstFilledRow: Integer;
    function lastFilledRow: Integer;
    function firstFilledCol: Integer;
    function lastFilledCol: Integer;

    function getPrintFit(var wPages, hPages: Integer): boolean;
    procedure setPrintFit(wPages, hPages: Integer);

    procedure setPrintRepeatRows(rowFirst, rowLast: Integer);
    procedure setPrintRepeatCols(colFirst, colLast: Integer);
    procedure setPrintArea(rowFirst, rowLast, colFirst, colLast: Integer);
    
    procedure clearPrintRepeats();
    procedure clearPrintArea();

    function getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer): boolean; overload;
    function getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer; scopeId: Integer): boolean; overload;
    function getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer; scopeId: Integer; var hidden: boolean): boolean; overload;
    function setNamedRange(const name: PWideChar; rowFirst, rowLast, colFirst, colLast: Integer): boolean; overload;
    function setNamedRange(const name: PWideChar; rowFirst, rowLast, colFirst, colLast, scopeId: Integer): boolean; overload;
    function delNamedRange(const name: PWideChar): boolean; overload;
    function delNamedRange(const name: PWideChar; scopeId: Integer): boolean; overload;

    function namedRangeSize(): Integer;
    function namedRange(index: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): WideString;

    function tableSize(): Integer;
    function table(index: Integer; var rowFirst, rowLast, colFirst, colLast, headerRowCount, totalsRowCount: Integer): WideString;

    function hyperlinkSize(): Integer;
    function hyperlink(index: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): WideString;
    function delHyperlink(index: Integer): boolean;
    procedure addHyperlink(const hyperlink: PWideChar; rowFirst, rowLast, colFirst, colLast: Integer);

    function autoFilter(): TXLAutoFilter;
    procedure applyFilter();
    procedure removeFilter();

    procedure setProtect(protect: boolean; const password: PWideChar); overload;
    procedure setProtect(protect: boolean; const password: PWideChar; const prot: EnhancedProtection); overload;

    procedure getTopLeftView(var row, col: Integer);
    procedure setTopLeftView(row, col: Integer);
    procedure setAutoFitArea(rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer);

    procedure addrToRowCol(const addr: PWideChar; var row, col: Integer); overload;
    procedure addrToRowCol(const addr: PWideChar; var row, col: Integer; var rowRelative, colRelative: boolean); overload;
    function rowColToAddr(row, col: Integer): WideString; overload;
    function rowColToAddr(row, col: Integer; rowRelative, colRelative: boolean): WideString; overload;

    procedure setTabColor(color: Color); overload;
    procedure setTabColor(red: Integer; green: Integer; blue: Integer); overload;

    function addIgnoredError(rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer; iError: IgnoredError): boolean;

    procedure addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value: PWideChar); overload;
    procedure addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value1: PWideChar; const value2: PWideChar); overload;
    procedure addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value1: PWideChar; const value2: PWideChar;
                                allowBlank: boolean; hideDropDown: boolean; showInputMessage: boolean; showErrorMessage: boolean; const promptTitle: PWideChar; const prompt: PWideChar;
                                const errorTitle: PWideChar; const error: PWideChar; errorStyle: DataValidationErrorStyle); overload;

    procedure addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value: double); overload;
    procedure addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double); overload;
    procedure addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double;
                                allowBlank: boolean; hideDropDown: boolean; showInputMessage: boolean; showErrorMessage: boolean; const promptTitle: PWideChar; const prompt: PWideChar;
                                const errorTitle: PWideChar; const error: PWideChar; errorStyle: DataValidationErrorStyle); overload;

    procedure removeDataValidations();

  published
    constructor Create(handle: SheetHandle; book: BookHandle);
    destructor Destroy; override;
  end;

type
  TXLBook = class(TObject)
  private
    handle: BookHandle;

    fonts : array of TXLFont;
    formats: array of TXLFormat;
    richStrings: array of TXLRichString;
    sheets: array of TXLSheet;

    procedure Release;

    function GetRefR1C1: boolean;
    procedure SetRefR1C1(const value: boolean);

    function GetDate1904: boolean;
    procedure SetDate1904(const value: boolean);

    function GetTemplate: boolean;
    procedure SetTemplate(const value: boolean);

  protected
  public
    function load(filename: PWideChar): boolean;
    function save(filename: PWideChar): boolean;

    function loadRaw(data: PByteArray; size: Cardinal): boolean;
    function saveRaw(var Buffer: ByteArray; var size: Cardinal): boolean;

    function addSheet(name: PWideChar): TXLSheet; overload;
    function addSheet(name: PWideChar; initSheet: TXLSheet): TXLSheet; overload;
    function insertSheet(index: Integer; name: PWideChar): TXLSheet; overload;
    function insertSheet(index: Integer; name: PWideChar; initSheet: TXLSheet): TXLSheet; overload;
    function getSheet(index: Integer): TXLSheet;
    function sheetType(index: Integer): SheetType;
    function moveSheet(srcIndex: Integer; dstIndex: Integer): boolean;
    function delSheet(index: Integer): boolean;
    function sheetCount(): Integer;

    function addFormat(): TXLFormat; overload;
    function addFormat(initFormat: TXLFormat): TXLFormat; overload;
    function addFont(): TXLFont; overload;
    function addFont(initFont: TXLFont): TXLFont; overload;
    function addRichString(): TXLRichString;
    function addCustomNumFormat(customNumFormat: PWideChar): Integer;
    function customNumFormat(fmt: Integer): WideString;

    function format(index: Integer): TXLFormat;
    function formatSize(): Integer;

    function font(index: Integer): TXLFont;
    function fontSize(): Integer;

    function datePack(year: Integer; month: Integer; day: Integer): double; overload;
    function datePack(year: Integer; month: Integer; day: Integer; hour: Integer; min: Integer; sec: Integer): double; overload;
    function datePack(year: Integer; month: Integer; day: Integer; hour: Integer; min: Integer; sec: Integer; msec: Integer): double; overload;

    function dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer): boolean; overload;
    function dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer; var hour: Integer; var min: Integer; var sec: Integer): boolean; overload;
    function dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer; var hour: Integer; var min: Integer; var sec: Integer; var msec: Integer): boolean; overload;

    function colorPack(red: Integer; green: Integer; blue: Integer): Color;
    procedure colorUnpack(color: Color; var red: Integer; var green: Integer; var blue: Integer);

    function activeSheet(): Integer;
    procedure setActiveSheet(index: Integer);

    function pictureSize(): Integer;
    function getPicture(index: Integer; var Buffer: ByteArray; var size: Cardinal): PictureType;

    function addPicture(filename: PWideChar): Integer;
    function addPicture2(data: PByteArray; size: Cardinal): Integer;
    function addPictureAsLink(filename: PWideChar; insert: boolean): Integer;

    function defaultFont(var fontSize: Integer): WideString;
    procedure setDefaultFont(fontName: PWideChar; fontSize: Integer);

    property refR1C1: boolean read GetRefR1C1 write SetRefR1C1;

    procedure setKey(name: PWideChar; key: PWideChar);

    function rgbMode(): boolean;
    procedure setRgbMode(rgbMode: boolean);

    function calcMode(): CalcModeType;
    procedure setCalcMode(calcMode: CalcModeType);

    function biffVersion(): Integer;

    property date1904: boolean read GetDate1904 write SetDate1904;
    property template: boolean read GetTemplate write SetTemplate;

    function setLocale(locale: PAnsiChar): boolean;

    function errorMessage(): WideString;

  published
    destructor Destroy; override;
  end;


  TBinBook = class(TXLBook)
  published
    constructor Create;    
  end;

  TXmlBook = class(TXLBook)
  published
    constructor Create;    
  end;
  

implementation

constructor TXLRichString.Create(handle: RichStringHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLRichString.Destroy;
var
  i : Integer;
begin
  for i := 0 to High(fonts) do fonts[i].Destroy;
  inherited Destroy;
end;

function TXLRichString.addFont(): TXLFont;
var
  fntHandle: FontHandle;
begin
  fntHandle := xlRichStringAddFont(handle, nil);
  if fntHandle <> nil then
  begin
    result := TXLFont.Create(fntHandle, handle);
    fonts := fonts + [result];
  end
  else
    result := nil;
end;

function TXLRichString.addFont(initFont: TXLFont): TXLFont;
var
  fntHandle: FontHandle;
begin
  fntHandle := xlRichStringAddFont(handle, initFont.handle);
  if fntHandle <> nil then
  begin
    result := TXLFont.Create(fntHandle, handle);
    fonts := fonts + [result];
  end
  else
    result := nil;
end;

procedure TXLRichString.addText(const text: PWideChar);
begin
   xlRichStringAddText(handle, text, nil);
end;

procedure TXLRichString.addText(const text: PWideChar; font: TXLFont);
begin
   xlRichStringAddText(handle, text, font.handle);
end;

function TXLRichString.getText(index: Integer): WideString;
var
   fntHandle: FontHandle;
begin
   result := WideString(xlRichStringGetText(handle, index, @fntHandle));
end;

function TXLRichString.getText(index: Integer; var font: TXLFont): WideString;
var
   fntHandle: FontHandle;
   ret: WideString;
begin
   ret := WideString(xlRichStringGetText(handle, index, @fntHandle));
   if fntHandle <> nil then
   begin
       font := TXLFont.Create(fntHandle, book);
       fonts := fonts + [font];
   end;
   result := ret;
end;

function TXLRichString.textSize(): Integer;
begin
  result := xlRichStringTextSize(handle);
end;

/////////

constructor TXLFilterColumn.Create(handle: FilterColumnHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLFilterColumn.Destroy;
begin
  inherited Destroy;
end;

function TXLFilterColumn.index(): Integer;
begin
  result := xlFilterColumnIndex(handle);
end;

function TXLFilterColumn.filterType(): Filter;
begin
  result := LibXL.Filter(xlFilterColumnFilterType(handle));
end;

function TXLFilterColumn.filterSize(): Integer;
begin
  result := xlFilterColumnFilterSize(handle);
end;

function TXLFilterColumn.filter(index: Integer): PWideChar;
begin
  result := xlFilterColumnFilter(handle, index);
end;

procedure TXLFilterColumn.addFilter(const value: PWideChar);
begin
  xlFilterColumnAddFilter(handle, value);
end;

function TXLFilterColumn.getTop10(var value: double; var top, percent: boolean): boolean;
var
  iTop, iPercent: Integer;
begin
  result := xlFilterColumnGetTop10(handle, @value, @iTop, @iPercent) > 0;
  top := iTop > 0;
  percent := iPercent > 0;
end;

procedure TXLFilterColumn.setTop10(value: double);
begin
  xlFilterColumnSetTop10(handle, value, 1, 0);
end;

procedure TXLFilterColumn.setTop10(value: double; top: boolean);
begin
  xlFilterColumnSetTop10(handle, value, Integer(top), 0);
end;

procedure TXLFilterColumn.setTop10(value: double; top, percent: boolean);
begin
  xlFilterColumnSetTop10(handle, value, Integer(top), Integer(percent));
end;

function TXLFilterColumn.getCustomFilter(var op1: Operator; var v1: PWideChar; var op2: Operator; var v2: PWideChar; var andOp: boolean): boolean;
var
  iOp1, iOp2, iAndOp: Integer;
begin
  result := xlFilterColumnGetCustomFilter(handle, @iOp1, @v1, @iOp2, @v2, @iAndOp) > 0;
  op1 := Operator(iOp1);
  op2 := Operator(iOp2);
  andOp := iAndOp > 0;
end;

procedure TXLFilterColumn.setCustomFilter(op: Operator; const v: PWideChar);
begin
  xlFilterColumnSetCustomFilter(handle, Integer(op), v);
end;

procedure TXLFilterColumn.setCustomFilter(op1: Operator; const v1: PWideChar; op2: Operator; const v2: PWideChar; andOp: boolean);
begin
  xlFilterColumnSetCustomFilterEx(handle, Integer(op1), v1, Integer(op2), v2, Integer(andOp));
end;

procedure TXLFilterColumn.clear();
begin
  xlFilterColumnClear(handle);
end;

constructor TXLAutoFilter.Create(handle: AutoFilterHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLAutoFilter.Destroy;
var
  i : Integer;
begin
  for i := 0 to High(filterColumns) do filterColumns[i].Destroy;
  inherited Destroy;
end;

function TXLAutoFilter.getRef(var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
begin
  result := xlAutoFilterGetRef(handle, @rowFirst, @rowLast, @colFirst, @colLast) > 0;
end;

procedure TXLAutoFilter.setRef(rowFirst, rowLast, colFirst, colLast: Integer);
begin
  xlAutoFilterSetRef(handle, rowFirst, rowLast, colFirst, colLast);
end;

function TXLAutoFilter.column(colId: Integer): TXLFilterColumn;
var
  fcHandle: FilterColumnHandle;
begin
  fcHandle := xlAutoFilterColumn(handle, colId);
  result := TXLFilterColumn.Create(fcHandle, book);
  filterColumns := filterColumns + [result];
end;

function TXLAutoFilter.columnSize(): Integer;
begin
  result := xlAutoFilterColumnSize(handle);
end;

function TXLAutoFilter.columnByIndex(index: Integer): TXLFilterColumn;
var
  fcHandle:  FilterColumnHandle;
begin
  fcHandle := xlAutoFilterColumnByIndex(handle, index);
  result := TXLFilterColumn.Create(fcHandle, book);
  filterColumns := filterColumns + [result];
end;

function TXLAutoFilter.getSortRange(var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
begin
  result := xlAutoFilterGetSortRange(handle, @rowFirst, @rowLast, @colFirst, @colLast) > 0;
end;

function TXLAutoFilter.getSort(var columnIndex: Integer; var descending: boolean): boolean;
var
  iDescending: Integer;
begin
  result := xlAutoFilterGetSort(handle, @columnIndex, @iDescending) > 0;
  descending := iDescending > 0;
end;

function TXLAutoFilter.setSort(columnIndex: Integer; descending: boolean): boolean;
begin
  result := xlAutoFilterSetSort(handle, columnIndex, Integer(descending)) > 0;
end;

constructor TXLFont.Create(handle: FontHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLFont.Destroy;
begin
  inherited Destroy;
end;

function TXLFont.Getbold: boolean;
begin
  result := (xlFontBold(handle) > 0);
end;

procedure TXLFont.Setbold(const value: boolean);
begin
  xlFontSetBold(handle, Integer(value));
end;

function TXLFont.GetName: WideString;
begin
  result := WideString(xlFontName(handle));
end;

procedure TXLFont.SetName(const value: WideString);
begin
  xlFontSetName(handle, @(value[1]));  
end;

function TXLFont.Getsize: Integer;
begin
  result := xlFontSize(handle);
end;

procedure TXLFont.Setsize(const value: Integer);
begin
  xlFontSetSize(handle, value);
end;

function TXLFont.GetUnderline: Underline;
begin
  result := LibXL.Underline(xlFontUnderline(handle));
end;

procedure TXLFont.SetUnderline(const value: LibXL.Underline);
begin
  xlFontSetUnderline(handle, Integer(value));
end;

function TXLFont.Getcolor: Color;
begin
  result := LibXL.Color(xlFontColor(handle));
end;

procedure TXLFont.Setcolor(const value: Color);
begin
  xlFontSetColor(handle, Integer(value));
end;

procedure TXLFont.Setitalic(const value: boolean);
begin
  xlFontSetItalic(handle, Integer(value));
end;

function TXLFont.Getitalic: boolean;
begin
  result := (xlFontItalic(handle) > 0);
end;

procedure TXLFont.Setscript(const value: Script);
begin
  xlFontSetScript(handle, Integer(value));
end;

function TXLFont.Getscript: Script;
begin
  result := LibXL.Script(xlFontScript(handle));
end;

procedure TXLFont.SetstrikeOut(const value: boolean);
begin
  xlFontSetStrikeOut(handle, Integer(value));
end;

function TXLFont.GetstrikeOut: boolean;
begin
  result := (xlFontStrikeOut(handle) > 0);
end;

constructor TXLFormat.Create(handle: FormatHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLFormat.Destroy;
begin
  if ffont <> nil then ffont.Destroy;
  inherited Destroy;
end;

procedure TXLFormat.SetborderDiagonal(const value: BorderDiagonal);
begin
  xlFormatSetBorderDiagonal(handle, Integer(value));
end;

procedure TXLFormat.SetBorderDiagonalStyle(const value: BorderStyle);
begin
  xlFormatSetBorderDiagonalStyle(handle, Integer(value));
end;

procedure TXLFormat.SetborderDiagonalColor(const value: Color);
begin
  xlFormatSetBorderDiagonalColor(handle, Integer(value));
end;

procedure TXLFormat.SetborderLeft(const value: BorderStyle);
begin
  xlFormatSetBorderLeft(handle, Integer(value));
end;

procedure TXLFormat.SetborderLeftColor(const value: Color);
begin
  xlFormatSetBorderLeftColor(handle, Integer(value));
end;

procedure TXLFormat.SetborderRight(const value: BorderStyle);
begin
  xlFormatSetBorderRight(handle, Integer(value));
end;

procedure TXLFormat.SetborderRightColor(const value: Color);
begin
  xlFormatSetBorderRightColor(handle, Integer(value));
end;

procedure TXLFormat.SetborderTop(const value: BorderStyle);
begin
  xlFormatSetBorderTop(handle, Integer(value));
end;

procedure TXLFormat.SetborderTopColor(const value: Color);
begin
  xlFormatSetBorderTopColor(handle, Integer(value));
end;

procedure TXLFormat.SetfillPattern(const value: FillPattern);
begin
  xlFormatSetFillPattern(handle, Integer(value));
end;

procedure TXLFormat.Setfont(const value: TXLFont);
begin
  xlFormatSetFont(handle, value.handle);  
end;

function TXLFormat.GetalignH: AlignH;
begin
  result := LibXL.AlignH(xlFormatAlignH(handle));
end;

function TXLFormat.GetAlignV: AlignV;
begin
  result := LibXL.AlignV(xlFormatAlignV(handle));
end;

function TXLFormat.GetborderBottom: BorderStyle;
begin
  result := BorderStyle(xlFormatBorderBottom(handle));
end;

function TXLFormat.GetborderBottomColor: Color;
begin
  result := Color(xlFormatBorderBottomColor(handle));
end;

function TXLFormat.GetborderDiagonal: BorderDiagonal;
begin
  result := LibXL.BorderDiagonal(xlFormatBorderDiagonal(handle));
end;

function TXLFormat.GetBorderDiagonalStyle: BorderStyle;
begin
  result := LibXL.BorderStyle(xlFormatBorderDiagonalStyle(handle));
end;

function TXLFormat.GetborderDiagonalColor: Color;
begin
  result := Color(xlFormatBorderDiagonalColor(handle));
end;

function TXLFormat.GetborderLeft: BorderStyle;
begin
  result := BorderStyle(xlFormatBorderLeft(handle));
end;

function TXLFormat.GetborderLeftColor: Color;
begin
  result := Color(xlFormatBorderLeftColor(handle));
end;

function TXLFormat.GetborderRight: BorderStyle;
begin
  result := BorderStyle(xlFormatBorderRight(handle));
end;

function TXLFormat.GetborderRightColor: Color;
begin
  result := Color(xlFormatBorderRightColor(handle));
end;

function TXLFormat.GetborderTop: BorderStyle;
begin
  result := BorderStyle(xlFormatBorderTop(handle));
end;

function TXLFormat.GetborderTopColor: Color;
begin
  result := Color(xlFormatBorderTopColor(handle));
end;

function TXLFormat.GetfillPattern: FillPattern;
begin
  result := LibXL.FillPattern(xlFormatFillPattern(handle));
end;

procedure TXLFormat.Sethidden(const value: boolean);
begin
  xlFormatSetHidden(handle, Integer(value));
end;

procedure TXLFormat.Setindent(const value: Integer);
begin
  xlFormatSetIndent(handle, Integer(value));
end;

procedure TXLFormat.Setlocked(const value: boolean);
begin
  xlFormatSetLocked(handle, Integer(value));
end;

function TXLFormat.Getfont: TXLFont;
var
  f: FontHandle;
begin
  f := xlFormatFont(handle);
  if f <> nil then
  begin
    result := TXLFont.Create(f, book);
    ffont := result;
  end
  else
    result := nil;
end;

function TXLFormat.Gethidden: boolean;
begin
  result := (xlFormatHidden(handle) > 0);
end;

function TXLFormat.Getindent: Integer;
begin
  result := (xlFormatIndent(handle));
end;

function TXLFormat.Getlocked: boolean;
begin
  result := (xlFormatLocked(handle) > 0);
end;

function TXLFormat.GetpatternBackgroundColor: Color;
begin
  result := Color(xlFormatPatternBackgroundColor(handle));
end;

function TXLFormat.GetpatternForegroundColor: Color;
begin
  result := Color(xlFormatPatternForegroundColor(handle));
end;

function TXLFormat.Getrotation: Integer;
begin
  result := xlFormatRotation(handle);
end;

function TXLFormat.GetshrinkToFit: boolean;
begin
  result := (xlFormatShrinkToFit(handle) > 0);
end;

function TXLFormat.Getwrap: boolean;
begin
  result := (xlFormatWrap(handle) > 0);
end;

function TXLFormat.NumFormat: Integer;
begin
  result := xlFormatNumFormat(handle);
end;

procedure TXLFormat.SetalignH(const value: AlignH);
begin
  xlFormatSetAlignH(handle, Integer(value));
end;

procedure TXLFormat.SetAlignV(const value: AlignV);
begin
  xlFormatSetAlignV(handle, Integer(value));
end;

procedure TXLFormat.setNumFormat(NumFormat: Integer);
begin
  xlFormatSetNumFormat(handle, NumFormat);
end;

procedure TXLFormat.setNumFormat(value: NumFormat);
begin
  xlFormatSetNumFormat(handle, Integer(value));
end;

procedure TXLFormat.SetpatternBackgroundColor(const value: Color);
begin
  xlFormatSetPatternBackgroundColor(handle, Integer(value));
end;

procedure TXLFormat.SetpatternForegroundColor(const value: Color);
begin
  xlFormatSetPatternForegroundColor(handle, Integer(value));
end;

procedure TXLFormat.Setrotation(const value: Integer);
begin
  xlFormatSetRotation(handle, value);
end;

procedure TXLFormat.SetshrinkToFit(const value: boolean);
begin
  xlFormatSetShrinkToFit(handle, Integer(value));
end;

procedure TXLFormat.Setwrap(const value: boolean);
begin
  xlFormatSetWrap(handle, Integer(value));
end;

procedure TXLFormat.setBorder(style: BorderStyle);
begin
  xlFormatSetBorder(handle, Integer(style));
end;

procedure TXLFormat.SetborderBottom(const value: BorderStyle);
begin
  xlFormatSetBorderBottom(handle, Integer(value));
end;

procedure TXLFormat.SetborderBottomColor(const value: Color);
begin
  xlFormatSetBorderBottomColor(handle, Integer(value));
end;

procedure TXLFormat.setBorderColor(value: Color);
begin
  xlFormatSetBorderColor(handle, Integer(value));
end;


constructor TXLSheet.Create(handle: SheetHandle; book: BookHandle);
begin
  inherited Create;
  self.book := book;
  self.handle := handle;
end;

destructor TXLSheet.Destroy;
var
  i: Integer;
begin
  for i := 0 to High(richStrings) do richStrings[i].Destroy;
  for i := 0 to High(formats) do formats[i].Destroy;
  if xlAutoFilter <> nil then xlAutoFilter.Destroy;
  inherited Destroy;
end;

function TXLSheet.ErrorMessage: WideString;
begin
  result := WideString(xlBookErrorMessage(book));
end;

function TXLSheet.getCellType(row, col: Integer): CellType;
begin
  result := CellType(xlSheetCellType(handle, row, col));
end;

function TXLSheet.cellFormat(row, col: Integer): TXLFormat;
begin
  result := TXLFormat.Create(xlSheetCellFormat(handle, row, col), book);
  formats := formats + [result];
end;

procedure TXLSheet.setCellFormat(row: Integer; col: Integer; format: TXLFormat);
begin
  xlSheetSetCellFormat(handle, row, col, format.handle);
end;

function TXLSheet.isFormula(row, col: Integer): boolean;
var
  res: Integer;
begin
  res := xlSheetIsFormula(handle, row, col);
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.readStr(row, col: Integer): WideString;
var
  format: FormatHandle;
  s: PWideChar;
begin
  s := xlSheetReadStr(handle, row, col, @format);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

function TXLSheet.readStr(row, col: Integer; var format: TXLFormat): WideString;
var
  fmtHandle: FormatHandle;
  s: PWideChar;
begin
  s := xlSheetReadStr(handle, row, col, @fmtHandle);
  if s <> nil then
  begin
    if fmtHandle <> nil then
    begin
      format := TXLFormat.Create(fmtHandle, book);
      formats := formats + [format];
    end
    else
      format := nil;
    result := WideString(s)
  end
  else
    result := '';
end;

function TXLSheet.writeStr(row, col: Integer; value: PWideChar): boolean;
begin
  result := xlSheetWriteStr(handle, row, col, value, nil) > 0;
end;

function TXLSheet.writeStr(row, col: Integer; value: PWideChar; format: TXLFormat): boolean;
begin
  result := xlSheetWriteStr(handle, row, col, value, format.handle) > 0;
end;

function TXLSheet.readRichStr(row: Integer; col: Integer): TXLRichString;
var
  richString: RichStringHandle;
  format: FormatHandle;
begin
  richString := xlSheetReadRichStr(handle, row, col, @format);
  if richString <> nil then
  begin
    result := TXLRichString.Create(richString, handle);
    richStrings := richStrings + [result];
  end
  else
    result := nil;
end;

function TXLSheet.readRichStr(row: Integer; col: Integer; var format: TXLFormat): TXLRichString;
var
  richString: RichStringHandle;
begin
  richString := xlSheetReadRichStr(handle, row, col, @format.handle);
  if richString <> nil then
  begin
    result := TXLRichString.Create(richString, handle);
    richStrings := richStrings + [result];
  end
  else
    result := nil;
end;

function TXLSheet.writeRichStr(row: Integer; col: Integer; value: TXLRichString): boolean;
begin
  result := xlSheetWriteRichStr(handle, row, col, value.handle, nil) > 0;
end;

function TXLSheet.writeRichStr(row: Integer; col: Integer; value: TXLRichString; format: TXLFormat): boolean;
begin
  result := xlSheetWriteRichStr(handle, row, col, value.handle, format.handle) > 0;
end;

function TXLSheet.readNum(row, col: Integer): double;
var
  format: FormatHandle;
begin
  result := xlSheetReadNum(handle, row, col, @format);
end;

function TXLSheet.readNum(row, col: Integer; var format: TXLFormat): double;
var
  fmtHandle: FormatHandle;
begin
  result := xlSheetReadNum(handle, row, col, @fmtHandle);
  if fmtHandle <> nil then
  begin
    format := TXLFormat.Create(fmtHandle, book);
    formats := formats + [format];
  end
  else
    format := nil;
end;

function TXLSheet.writeNum(row, col: Integer; value: double): boolean;
begin
  result := xlSheetWriteNum(handle, row, col, value, nil) > 0;
end;

function TXLSheet.writeNum(row, col: Integer; value: double; format: TXLFormat): boolean;
begin
  result := xlSheetWriteNum(handle, row, col, value, format.handle) > 0;
end;

function TXLSheet.readBool(row, col: Integer): boolean;
var
  format: FormatHandle;
  res: Integer;
begin
  res := xlSheetReadBool(handle, row, col, @format);
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.readBool(row, col: Integer; var format: TXLFormat): boolean;
var
  fmtHandle: FormatHandle;
  res: Integer;
begin
  res := xlSheetReadBool(handle, row, col, @fmtHandle);
  if fmtHandle <> nil then
  begin
    format := TXLFormat.Create(fmtHandle, book);
    formats := formats + [format];
  end
  else
    format := nil;
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.writeBool(row, col: Integer; value: boolean): boolean;
begin
  result := xlSheetWriteBool(handle, row, col, Integer(value), nil) > 0;
end;

function TXLSheet.writeBool(row, col: Integer; value: boolean; format: TXLFormat): boolean;
begin
  result := xlSheetWriteBool(handle, row, col, Integer(value), format.handle) > 0;
end;

function TXLSheet.readBlank(row, col: Integer; var format: TXLFormat): boolean;
var
  fmtHandle: FormatHandle;
begin
  result := xlSheetReadBlank(handle, row, col, @fmtHandle) > 0;
  if fmtHandle <> nil then
  begin
    format := TXLFormat.Create(fmtHandle, book);
    formats := formats + [format];
  end
  else
    format := nil;  
end;

function TXLSheet.writeBlank(row, col: Integer; format: TXLFormat): boolean;
begin
  result := xlSheetWriteBlank(handle, row, col, format.handle) > 0;
end;

function TXLSheet.readFormula(row, col: Integer): WideString;
var
  format: FormatHandle;
  s: PWideChar;
begin
  s := xlSheetReadFormula(handle, row, col, @format);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

function TXLSheet.readFormula(row, col: Integer; var format: TXLFormat): WideString;
var
  fmtHandle: FormatHandle;
  s: PWideChar;
begin
  s := xlSheetReadFormula(handle, row, col, @fmtHandle);
  if s <> nil then
  begin
    if fmtHandle <> nil then
    begin
      format := TXLFormat.Create(fmtHandle, book);
      formats := formats + [format];
    end
    else
      format := nil;
    result := WideString(s);
  end
  else
    result := '';
end;

function TXLSheet.writeFormula(row, col: Integer; value: PWideChar): boolean;
begin
  result := xlSheetWriteFormula(handle, row, col, value, nil) > 0;
end;

function TXLSheet.writeFormula(row, col: Integer; value: PWideChar; format: TXLFormat): boolean;
begin
  result := xlSheetWriteFormula(handle, row, col, value, format.handle) > 0;
end;

function TXLSheet.writeFormulaNum(row, col: Integer; expr: PWideChar; value: double): boolean;
begin
  result := xlSheetWriteFormulaNum(handle, row, col, expr, value, nil) > 0;
end;

function TXLSheet.writeFormulaNum(row, col: Integer; expr: PWideChar; value: double; format: TXLFormat): boolean;
begin
  result := xlSheetWriteFormulaNum(handle, row, col, expr, value, format.handle) > 0;
end;

function TXLSheet.writeFormulaStr(row, col: Integer; expr: PWideChar; value: PWideChar): boolean;
begin
  result := xlSheetWriteFormulaStr(handle, row, col, expr, value, nil) > 0;
end;

function TXLSheet.writeFormulaStr(row, col: Integer; expr: PWideChar; value: PWideChar; format: TXLFormat): boolean;
begin
  result := xlSheetWriteFormulaStr(handle, row, col, expr, value, format.handle) > 0;
end;

function TXLSheet.writeFormulaBool(row, col: Integer; expr: PWideChar; value: Integer): boolean;
begin
  result := xlSheetWriteFormulaBool(handle, row, col, expr, value, nil) > 0;
end;

function TXLSheet.writeFormulaBool(row, col: Integer; expr: PWideChar; value: Integer; format: TXLFormat): boolean;
begin
  result := xlSheetWriteFormulaBool(handle, row, col, expr, value, format.handle) > 0;
end;

function TXLSheet.readComment(row, col: Integer): WideString;
var
  s: PWideChar;
begin
  s := xlSheetReadComment(handle, row, col);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLSheet.writeComment(row, col: Integer; value, author: PWideChar; width, height: Integer);
begin
  xlSheetWriteComment(handle, row, col, value, author, width, height);
end;

procedure TXLSheet.removeComment(row, col: Integer);
begin
  xlSheetRemoveComment(handle, row, col);
end;

function TXLSheet.isDate(row, col: Integer): boolean;
begin
  if xlSheetIsDate(handle, row, col) > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.isRichStr(row, col: Integer): boolean;
begin
  if xlSheetIsRichStr(handle, row, col) > 0 then
    result := true
  else
    result := false;
end;


function TXLSheet.readError(row, col: Integer): ErrorType;
begin
  result := ErrorType(xlSheetReadError(handle, row, col));
end;

procedure TXLSheet.writeError(row, col: Integer; error: ErrorType);
begin
  xlSheetWriteError(handle, row, col, Integer(error), nil);
end;

procedure TXLSheet.writeError(row, col: Integer; error: ErrorType; format: TXLFormat);
begin
  xlSheetWriteError(handle, row, col, Integer(error), format.handle);
end;

function TXLSheet.colWidth(col: Integer): double;
begin
  result := xlSheetColWidth(handle, col);
end;

function TXLSheet.rowHeight(col: Integer): double;
begin
  result := xlSheetRowHeight(handle, col);
end;

function TXLSheet.setCol(col: Integer; width: double; format: TXLFormat): boolean;
begin
  result := setCol(col, col, width, format);
end;

function TXLSheet.setCol(colFirst, colLast: Integer; width: double): boolean;
begin
  result := xlSheetSetCol(handle, colFirst, colLast, width, nil, 0) > 0;
end;

function TXLSheet.setCol(col: Integer; width: double): boolean;
begin
  result := setCol(col, col, width);
end;

function TXLSheet.setCol(colFirst, colLast: Integer; width: double; format: TXLFormat; hidden: boolean): boolean;
var
  fmtHandle: FormatHandle;
begin
  if format <> nil then
    fmtHandle := format.handle
  else
    fmtHandle := nil;
  result := xlSheetSetCol(handle, colFirst, colLast, width, fmtHandle, Integer(hidden)) > 0;
end;

function TXLSheet.setCol(col: Integer; width: double; format: TXLFormat; hidden: boolean): boolean;
begin
  result := setCol(col, col, width, format, hidden);
end;

function TXLSheet.setCol(colFirst, colLast: Integer; width: double; format: TXLFormat): boolean;
var
  fmtHandle: FormatHandle;
begin
  if format <> nil then
    fmtHandle := format.handle
  else
    fmtHandle := nil;
  result := xlSheetSetCol(handle, colFirst, colLast, width, fmtHandle, 0) > 0;
end;

function TXLSheet.setRow(row: Integer; height: double): boolean;
begin
  result := xlSheetSetRow(handle, row, height, nil, 0) > 0;
end;

function TXLSheet.setRow(row: Integer; height: double; format: TXLFormat): boolean;
begin
  result := xlSheetSetRow(handle, row, height, format.handle, 0) > 0;
end;

function TXLSheet.setRow(row: Integer; height: double; format: TXLFormat; hidden: boolean): boolean;
var
  fmtHandle: FormatHandle;
begin
  if format <> nil then
    fmtHandle := format.handle
  else
    fmtHandle := nil;
  result := xlSheetSetRow(handle, row, height, fmtHandle, Integer(hidden)) > 0;
end;

function TXLSheet.rowHidden(row: Integer): boolean;
begin
  result := xlSheetRowHidden(handle, row) > 0;
end;

function TXLSheet.setRowHidden(row: Integer; hidden: boolean): boolean;
begin
  result := xlSheetSetRowHidden(handle, row, Integer(hidden)) > 0;
end;

function TXLSheet.colHidden(col: Integer): boolean;
begin
  result := xlSheetColHidden(handle, col) > 0;
end;

function TXLSheet.setColHidden(col: Integer; hidden: boolean): boolean;
begin
  result := xlSheetSetColHidden(handle, col, Integer(hidden)) > 0;
end;

function TXLSheet.delMerge(row, col: Integer): boolean;
var
  res: Integer;
begin
  res := xlSheetDelMerge(handle, row, col);
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.getMerge(row, col: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
var
  res: Integer;
begin
  res := xlSheetGetMerge(handle, row, col, @rowFirst, @rowLast, @colFirst, @colLast);
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.setMerge(rowFirst, rowLast, colFirst, colLast: Integer): boolean;
var
  res: Integer;
begin
  res := xlSheetSetMerge(handle, rowFirst, rowLast, colFirst, colLast);
  if res > 0 then
    result := true
  else
    result := false;
end;

function TXLSheet.pictureSize: Integer;
begin
  result := xlSheetPictureSize(handle);
end;

function TXLSheet.getPicture(index: Integer): Integer;
var
  rowTop, colLeft, rowBottom, colRight, width, height, offset_x, offset_y: Integer;
begin
  result := xlSheetGetPicture(handle, index, @rowTop, @colLeft, @rowBottom, @colRight, @width, @height, @offset_x, @offset_y);
end;

function TXLSheet.getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight: Integer): Integer;
var
  width, height, offset_x, offset_y: Integer;
begin
  result := xlSheetGetPicture(handle, index, @rowTop, @colLeft, @rowBottom, @colRight, @width, @height, @offset_x, @offset_y);
end;

function TXLSheet.getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight, width, height: Integer): Integer;
var
  offset_x, offset_y: Integer;
begin
  result := xlSheetGetPicture(handle, index, @rowTop, @colLeft, @rowBottom, @colRight, @width, @height, @offset_x, @offset_y);
end;

function TXLSheet.getPicture(index: Integer; var rowTop, colLeft, rowBottom, colRight, width, height, offset_x, offset_y: Integer): Integer;
begin
  result := xlSheetGetPicture(handle, index, @rowTop, @colLeft, @rowBottom, @colRight, @width, @height, @offset_x, @offset_y);
end;

function TXLSheet.delHorPageBreak(row: Integer): boolean;
begin
  result := xlSheetSetHorPageBreak(handle, row, 0) > 0;
end;

function TXLSheet.delVerPageBreak(col: Integer): boolean;
begin
  result := xlSheetSetVerPageBreak(handle, col, 0) > 0;
end;

function TXLSheet.setHorPageBreak(row: Integer): boolean;
begin
  result := xlSheetSetHorPageBreak(handle, row, 1) > 0;
end;

function TXLSheet.setVerPageBreak(col: Integer): boolean;
begin
  result := xlSheetSetVerPageBreak(handle, col, 1) > 0;
end;

procedure TXLSheet.setPicture(row, col, pictureId, width, height: Integer);
begin
  xlSheetSetPicture2(handle, row, col, pictureId, width, height, 0, 0, 0);
end;

procedure TXLSheet.setPicture(row, col, pictureId: Integer; scale: double);
begin
  xlSheetSetPicture(handle, row, col, pictureId, scale, 0, 0, 0);
end;

procedure TXLSheet.setPicture(row, col, pictureId: Integer);
begin
  setPicture(row, col, pictureId, 1);
end;

function TXLSheet.removePictureByIndex(index: Integer): boolean;
begin
  result := xlSheetRemovePictureByIndex(handle, index) > 0;
end;

function TXLSheet.removePicture(row, col: Integer): boolean;
begin
  result := xlSheetRemovePicture(handle, row, col) > 0;
end;

function TXLSheet.getHorPageBreak(index: Integer): Integer;
begin
  result := xlSheetGetHorPageBreak(handle, index);
end;

function TXLSheet.getHorPageBreakSize: Integer;
begin
  result := xlSheetGetHorPageBreakSize(handle);
end;

function TXLSheet.getVerPageBreak(index: Integer): Integer;
begin
  result := xlSheetGetVerPageBreak(handle, index);
end;

function TXLSheet.getVerPageBreakSize: Integer;
begin
  result := xlSheetGetVerPageBreakSize(handle);
end;

procedure TXLSheet.split(row, col: Integer);
begin
  xlSheetSplit(handle, row, col);
end;

function TXLSheet.groupCols(colFirst, colLast: Integer; collapsed: boolean): boolean;
begin
  result := xlSheetGroupCols(handle, colFirst, colLast, Integer(collapsed)) > 0;
end;

function TXLSheet.groupRows(rowFirst, rowLast: Integer; collapsed: boolean): boolean;
begin
  result := xlSheetGroupRows(handle, rowFirst, rowLast, Integer(collapsed)) > 0;
end;

function TXLSheet.clear(rowFirst, rowLast, colFirst, colLast: Integer): boolean;
begin
  result := xlSheetClear(handle, rowFirst, rowLast, colFirst, colLast) > 0;
end;

function TXLSheet.insertCol(colFirst, colLast: Integer): boolean;
begin
  result := xlSheetInsertCol(handle, colFirst, colLast) > 0;
end;

function TXLSheet.insertRow(rowFirst, rowLast: Integer): boolean;
begin
  result := xlSheetInsertRow(handle, rowFirst, rowLast) > 0;
end;

function TXLSheet.removeCol(colFirst, colLast: Integer): boolean;
begin
  result := xlSheetRemoveCol(handle, colFirst, colLast) > 0;
end;

function TXLSheet.removeRow(rowFirst, rowLast: Integer): boolean;
begin
  result := xlSheetRemoveRow(handle, rowFirst, rowLast) > 0;
end;

function TXLSheet.copyCell(rowSrc, colSrc, rowDst, colDst: Integer): boolean;
begin
  result := xlSheetCopyCell(handle, rowSrc, colSrc, rowDst, colDst) > 0;
end;

procedure TXLSheet.setPrintRepeatRows(rowFirst, rowLast: Integer);
begin
  xlSheetSetPrintRepeatRows(handle, rowFirst, rowLast);
end;

procedure TXLSheet.setPrintRepeatCols(colFirst, colLast: Integer);
begin
  xlSheetSetPrintRepeatCols(handle, colFirst, colLast);
end;

procedure TXLSheet.setPrintArea(rowFirst, rowLast, colFirst, colLast: Integer);
begin
  xlSheetSetPrintArea(handle, rowFirst, rowLast, colFirst, colLast);
end;
    
procedure TXLSheet.clearPrintRepeats();
begin
  xlSheetClearPrintRepeats(handle);
end;

procedure TXLSheet.clearPrintArea();
begin
  xlSheetClearPrintArea(handle);
end;

function TXLSheet.getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer): boolean;
var
  hidden: boolean;
begin
  result := xlSheetGetNamedRange(handle, name, @rowFirst, @rowLast, @colFirst, @colLast, -2, @hidden) > 0;
end;

function TXLSheet.getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer; scopeId: Integer): boolean;
var
  hidden: boolean;
begin
  result := xlSheetGetNamedRange(handle, name, @rowFirst, @rowLast, @colFirst, @colLast, scopeId, @hidden) > 0;
end;

function TXLSheet.getNamedRange(const name: PWideChar; var rowFirst, rowLast, colFirst, colLast: Integer; scopeId: Integer; var hidden: boolean): boolean;
begin
  result := xlSheetGetNamedRange(handle, name, @rowFirst, @rowLast, @colFirst, @colLast, scopeId, @hidden) > 0;
end;

function TXLSheet.setNamedRange(const name: PWideChar; rowFirst, rowLast, colFirst, colLast: Integer): boolean;
begin
  result := xlSheetSetNamedRange(handle, name, rowFirst, rowLast, colFirst, colLast, -2) > 0;
end;

function TXLSheet.setNamedRange(const name: PWideChar; rowFirst, rowLast, colFirst, colLast, scopeId: Integer): boolean;
begin
  result := xlSheetSetNamedRange(handle, name, rowFirst, rowLast, colFirst, colLast, scopeId) > 0;
end;

function TXLSheet.delNamedRange(const name: PWideChar): boolean;
begin
  result := xlSheetDelNamedRange(handle, name, -2) > 0;
end;

function TXLSheet.delNamedRange(const name: PWideChar; scopeId: Integer): boolean;
begin
  result := xlSheetDelNamedRange(handle, name, scopeId) > 0;
end;

function TXLSheet.namedRangeSize(): Integer;
begin
  result := xlSheetNamedRangeSize(handle);
end;

function TXLSheet.namedRange(index: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): WideString;
var
  scopeId: Integer;
  hidden: boolean;
begin
  result := xlSheetNamedRange(handle, index, @rowFirst, @rowLast, @colFirst, @colLast, @scopeId, @hidden);
end;

function TXLSheet.tableSize(): Integer;
begin
  result := xlSheetTableSize(handle);
end;

function TXLSheet.table(index: Integer; var rowFirst, rowLast, colFirst, colLast, headerRowCount, totalsRowCount: Integer): WideString;
begin
  result := xlSheetTable(handle, index, @rowFirst, @rowLast, @colFirst, @colLast, @headerRowCount, @totalsRowCount);
end;

function TXLSheet.hyperlinkSize(): Integer;
begin
  result := xlSheetHyperlinkSize(handle);
end;

function TXLSheet.hyperlink(index: Integer; var rowFirst, rowLast, colFirst, colLast: Integer): WideString;
begin
  result := xlSheetHyperlink(handle, index, @rowFirst, @rowLast, @colFirst, @colLast);
end;

function TXLSheet.delHyperlink(index: Integer): boolean;
begin
  result := xlSheetDelHyperlink(handle, index) > 0;
end;

procedure TXLSheet.addHyperlink(const hyperlink: PWideChar; rowFirst, rowLast, colFirst, colLast: Integer);
begin
  xlSheetAddHyperlink(handle, hyperlink, rowFirst, rowLast, colFirst, colLast);
end;

function TXLSheet.autoFilter(): TXLAutoFilter;
var
  autoFilter: AutoFilterHandle;
begin
  autoFilter := xlSheetAutoFilter(handle);
  if autoFilter <> nil then
  begin
    result := TXLAutoFilter.Create(xlSheetAutoFilter(handle), book);
    xlAutoFilter := result
  end
  else
    result := nil;
end;

procedure TXLSheet.applyFilter();
begin
  xlSheetApplyFilter(handle);
end;

procedure TXLSheet.removeFilter();
begin
  xlSheetRemoveFilter(handle);
end;

procedure TXLSheet.setProtect(protect: boolean; const password: PWideChar);
begin
  xlSheetSetProtectEx(handle, Integer(protect), password, -1);
end;

procedure TXLSheet.setProtect(protect: boolean; const password: PWideChar; const prot: EnhancedProtection);
begin
  xlSheetSetProtectEx(handle, Integer(protect), password, Integer(prot));
end;

function TXLSheet.firstRow: Integer;
begin
  result := xlSheetFirstRow(handle);
end;

function TXLSheet.lastRow: Integer;
begin
  result := xlSheetLastRow(handle);
end;

function TXLSheet.firstCol: Integer;
begin
  result := xlSheetFirstCol(handle);
end;

function TXLSheet.lastCol: Integer;
begin
  result := xlSheetLastCol(handle);
end;

function TXLSheet.firstFilledRow: Integer;
begin
  result := xlSheetFirstFilledRow(handle);
end;

function TXLSheet.lastFilledRow: Integer;
begin
  result := xlSheetLastFilledRow(handle);
end;

function TXLSheet.firstFilledCol: Integer;
begin
  result := xlSheetFirstFilledCol(handle);
end;

function TXLSheet.lastFilledCol: Integer;
begin
  result := xlSheetLastFilledCol(handle);
end;

function TXLSheet.getPrintFit(var wPages, hPages: Integer): boolean;
begin
  result := xlSheetGetPrintFit(handle, @wPages, @hPages) > 0;
end;

procedure TXLSheet.setPrintFit(wPages, hPages: Integer);
begin
  xlSheetSetPrintFit(handle, wPages, hPages);
end;

function TXLSheet.GetDisplayGridlines: boolean;
begin
  result := (xlSheetDisplayGridlines(handle) > 0);
end;

procedure TXLSheet.SetDisplayGridlines(const value: boolean);
begin
  xlSheetSetDisplayGridlines(handle, Integer(value));
end;

function TXLSheet.GetprintGridlines: boolean;
begin
  result := (xlSheetPrintGridlines(handle) > 0);
end;

procedure TXLSheet.SetprintGridlines(const value: boolean);
begin
  xlSheetSetPrintGridlines(handle, Integer(value));
end;

procedure TXLSheet.SetZoom(const value: Integer);
begin
  xlSheetSetZoom(handle, value);
end;

function TXLSheet.GetZoom: Integer;
begin
  result := xlSheetZoom(handle);
end;

function TXLSheet.GetprintZoom: Integer;
begin
  result := xlSheetPrintZoom(handle);
end;

procedure TXLSheet.SetprintZoom(const value: Integer);
begin
  xlSheetSetPrintZoom(handle, value);
end;

procedure TXLSheet.Setlandscape(const value: boolean);
begin
  xlSheetSetLandscape(handle, Integer(value));
end;

function TXLSheet.Getlandscape: boolean;
begin
  result := (xlSheetLandscape(handle) > 0);
end;

function TXLSheet.GetPaper: Paper;
begin
  result := LibXL.Paper(xlSheetPaper(handle));
end;

procedure TXLSheet.SetPaper(const value: Paper);
begin
  xlSheetSetPaper(handle, Integer(value));
end;

function TXLSheet.Getheader: WideString;
var
  s: PWideChar;
begin
  s := xlSheetHeader(handle);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLSheet.Setheader(const value: WideString);
var
  PValue: PWideChar;
begin
  PValue := @(value[1]);
  xlSheetSetHeader(handle, PValue, 0.5);
end;

function TXLSheet.GetheaderMargin: double;
begin
  result := xlSheetHeaderMargin(handle);
end;

procedure TXLSheet.SetheaderMargin(const value: double);
begin
  xlSheetSetHeader(handle, @(self.header[1]), value);
end;

function TXLSheet.Getfooter: WideString;
var
  s: PWideChar;
begin
  s := xlSheetFooter(handle);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLSheet.Setfooter(const value: WideString);
var
  PValue: PWideChar;
begin
  PValue := @(value[1]);
  xlSheetSetFooter(handle, PValue, 0.5);
end;

function TXLSheet.GetfooterMargin: double;
begin
  result := xlSheetHeaderMargin(handle);
end;

procedure TXLSheet.SetfooterMargin(const value: double);
begin
  xlSheetSetHeader(handle, @(self.footer[1]), value);
end;

procedure TXLSheet.SethCenter(const value: boolean);
begin
  xlSheetSetHCenter(handle, Integer(value));
end;

function TXLSheet.GethCenter: boolean;
begin
  result := (xlSheetHCenter(handle) > 0);
end;

procedure TXLSheet.SetvCenter(const value: boolean);
begin
  xlSheetSetVCenter(handle, Integer(value));
end;

function TXLSheet.GetvCenter: boolean;
begin
  result := (xlSheetVCenter(handle) > 0);
end;

procedure TXLSheet.SetmarginLeft(const value: double);
begin
  xlSheetSetMarginLeft(handle, value);
end;

function TXLSheet.GetmarginLeft: double;
begin
  result := xlSheetMarginLeft(handle);
end;

procedure TXLSheet.SetmarginRight(const value: double);
begin
  xlSheetSetMarginRight(handle, value);
end;

function TXLSheet.GetmarginRight: double;
begin
  result := xlSheetMarginRight(handle);
end;

procedure TXLSheet.SetmarginTop(const value: double);
begin
  xlSheetSetMarginTop(handle, value);
end;

function TXLSheet.GetmarginTop: double;
begin
  result := xlSheetMarginTop(handle);
end;

procedure TXLSheet.SetmarginBottom(const value: double);
begin
  xlSheetSetMarginBottom(handle, value);
end;

function TXLSheet.GetmarginBottom: double;
begin
  result := xlSheetMarginBottom(handle);
end;

function TXLSheet.Getname: WideString;
var
  s: PWideChar;
begin
  s := xlSheetName(handle);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLSheet.Setname(const value: WideString);
begin
  xlSheetSetName(handle, @(value[1]));
end;

function TXLSheet.GetprintRowCol: boolean;
begin
  result := (xlSheetPrintRowCol(handle) > 0);
end;

procedure TXLSheet.SetprintRowCol(const value: boolean);
begin
  xlSheetSetPrintRowCol(handle, Integer(value));
end;

function TXLSheet.GetprotectImpl: boolean;
begin
  result := (xlSheetProtect(handle) > 0);
end;

procedure TXLSheet.SetprotectImpl(const value: boolean);
begin
  xlSheetSetProtect(handle, Integer(value));
end;

function TXLSheet.Getrighttoleft: boolean;
begin
  result := (xlSheetRightToLeft(handle) > 0);
end;

procedure TXLSheet.Setrighttoleft(const value: boolean);
begin
  xlSheetSetRightToLeft(handle, Integer(value));
end;

function TXLSheet.GetHidden: SheetState;
begin
  result := SheetState(xlSheetHidden(handle));
end;

procedure TXLSheet.SetHidden(const value: SheetState);
begin
  xlSheetSetHidden(handle, Integer(value));
end;

procedure TXLSheet.getTopLeftView(var row, col: Integer);
begin
  xlSheetGetTopLeftView(handle, @row, @col);  
end;

procedure TXLSheet.setTopLeftView(row, col: Integer);
begin
  xlSheetSetTopLeftView(handle, row, col);  
end;

procedure TXLSheet.setAutoFitArea(rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer);
begin
  xlSheetSetAutoFitArea(handle, rowFirst, colFirst, rowLast, colLast);
end;

function TXLSheet.GetGroupSummaryBelow: boolean;
begin
  result := (xlSheetGroupSummaryBelow(handle) > 0);
end;

procedure TXLSheet.SetGroupSummaryBelow(const value: boolean);
begin
  xlSheetSetGroupSummaryBelow(handle, Integer(value));
end;

function TXLSheet.GetGroupSummaryRight: boolean;
begin
  result := (xlSheetGroupSummaryRight(handle) > 0);
end;

procedure TXLSheet.SetGroupSummaryRight(const value: boolean);
begin
  xlSheetSetGroupSummaryRight(handle, Integer(value));
end;

procedure TXLSheet.addrToRowCol(const addr: PWideChar; var row, col: Integer);
begin
  xlSheetAddrToRowCol(handle, addr, @row, @col, nil, nil);
end;

procedure TXLSheet.addrToRowCol(const addr: PWideChar; var row, col: Integer; var rowRelative, colRelative: boolean);
var
  localRowRelative, localColRelative: Integer;
begin
  xlSheetAddrToRowCol(handle, addr, @row, @col, @localRowRelative, @localColRelative);
  rowRelative := localRowRelative > 0;
  colRelative := localColRelative > 0;
end;

function TXLSheet.rowColToAddr(row, col: Integer): WideString;
var
  s: PWideChar;
begin
  s := xlSheetRowColToAddr(handle, row, col, 0, 0);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

function TXLSheet.rowColToAddr(row, col: Integer; rowRelative, colRelative: boolean): WideString;
var
  s: PWideChar;
begin
  s := xlSheetRowColToAddr(handle, row, col, Integer(rowRelative), Integer(colRelative));
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLSheet.setTabColor(color: Color);
begin
  xlSheetSetTabColor(handle, Integer(color));
end;

procedure TXLSheet.setTabColor(red: Integer; green: Integer; blue: Integer);
begin
  xlSheetSetTabRgbColor(handle, red, green, blue);
end;

function TXLSheet.addIgnoredError(rowFirst: Integer; colFirst: Integer; rowLast: Integer; colLast: Integer; iError: IgnoredError): boolean;
begin
  result := xlSheetAddIgnoredError(handle,  rowFirst,  colFirst, rowLast, colLast, Integer(iError)) > 0;
end;

procedure TXLSheet.addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value: PWideChar);
begin
  xlSheetAddDataValidation(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value, nil);
end;

procedure TXLSheet.addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value1: PWideChar; const value2: PWideChar);
begin
  xlSheetAddDataValidation(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value1, value2);
end;

procedure TXLSheet.addDataValidation(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; const value1: PWideChar; const value2: PWideChar;
                                allowBlank: boolean; hideDropDown: boolean; showInputMessage: boolean; showErrorMessage: boolean; const promptTitle: PWideChar; const prompt: PWideChar;
                                const errorTitle: PWideChar; const error: PWideChar; errorStyle: DataValidationErrorStyle);
begin
  xlSheetAddDataValidationEx(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value1, value2, Integer(allowBlank), Integer(hideDropDown), Integer(showInputMessage), Integer(showErrorMessage), promptTitle, prompt, errorTitle, error, Integer(errorStyle));
end;

procedure TXLSheet.addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value: double);
begin
  xlSheetAddDataValidationDouble(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value, 0);
end;

procedure TXLSheet.addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double);
begin
  xlSheetAddDataValidationDouble(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value1, value2);
end;

procedure TXLSheet.addDataValidationDouble(vtype: DataValidationType; op: DataValidationOperator; rowFirst: Integer; rowLast: Integer; colFirst: Integer; colLast: Integer; value1: double; value2: double;
                                allowBlank: boolean; hideDropDown: boolean; showInputMessage: boolean; showErrorMessage: boolean; const promptTitle: PWideChar; const prompt: PWideChar;
                                const errorTitle: PWideChar; const error: PWideChar; errorStyle: DataValidationErrorStyle);
begin
  xlSheetAddDataValidationDoubleEx(handle, Integer(vtype), Integer(op), rowFirst, rowLast, colFirst, colLast, value1, value2, Integer(allowBlank), Integer(hideDropDown), Integer(showInputMessage), Integer(showErrorMessage), promptTitle, prompt, errorTitle, error, Integer(errorStyle));
end;

procedure TXLSheet.removeDataValidations();
begin
  xlSheetRemoveDataValidations(handle);
end;

function TXLBook.load(filename: PWideChar): boolean;
begin
  result := xlBookLoad(handle, filename) > 0;
end;

function TXLBook.save(filename: PWideChar): boolean;
begin
  result := xlBookSave(handle, filename) > 0;
end;

function TXLBook.loadRaw(data: PByteArray; size: Cardinal): boolean;
begin
  result := xlBookLoadRaw(handle, addr(data^[0]), size) > 0;
end;

function TXLBook.saveRaw(var Buffer: ByteArray; var size: Cardinal): boolean;
var
  ptr: PByteArray;
  i: Integer;
  p: PByte;
begin
  if xlBookSaveRaw(handle, @ptr, @size) <> 0 then
  begin
    SetLength(Buffer, size);
    p := @(ptr^);
    for i := 0 to size - 1 do
    begin
      Buffer[i] := Byte(p^);
      Inc(p);
    end;
    result := true;
  end
  else
    result := false;
end;

function TXLBook.addSheet(name: PWideChar): TXLSheet;
var
  new_sheetHandle: SheetHandle;
begin
  new_sheetHandle := xlBookAddSheet(handle, name, nil);
  if new_sheetHandle <> nil then
  begin
    result := TXLSheet.Create(new_sheetHandle, self.handle);
    sheets := sheets + [result];
  end
  else
    result := nil;
end;

function TXLBook.addSheet(name: PWideChar; initSheet: TXLSheet): TXLSheet;
var
  new_sheetHandle: SheetHandle;
begin
  new_sheetHandle := xlBookAddSheet(handle, name, initSheet.handle);
  if new_sheetHandle <> nil then
  begin
    result := TXLSheet.Create(new_sheetHandle, self.handle);
    sheets := sheets + [result];
  end
  else
    result := nil;
end;

function TXLBook.insertSheet(index: Integer; name: PWideChar): TXLSheet;
var
  new_sheetHandle: SheetHandle;
begin
  new_sheetHandle := xlBookInsertSheet(handle, index, name, nil);
  if new_sheetHandle <> nil then
  begin
    result := TXLSheet.Create(new_sheetHandle, self.handle);
    sheets := sheets + [result];
  end
  else
    result := nil;
end;

function TXLBook.insertSheet(index: Integer; name: PWideChar; initSheet: TXLSheet): TXLSheet;
var
  new_sheetHandle: SheetHandle;
begin
  new_sheetHandle := xlBookInsertSheet(handle, index, name, initSheet);
  if new_sheetHandle <> nil then
  begin
    result := TXLSheet.Create(new_sheetHandle, self.handle);
    sheets := sheets + [result];
  end
  else 
    result := nil;
end;

function TXLBook.getSheet(index: Integer): TXLSheet;
var
  new_sheetHandle: SheetHandle;
begin
  new_sheetHandle := xlBookGetSheet(handle, index);
  if new_sheetHandle <> nil then
  begin
    result := TXLSheet.Create(new_sheetHandle, self.handle);
    sheets := sheets + [result];
  end
  else
    result := nil;
end;

function TXLBook.sheetType(index: Integer): SheetType;
var
  retval: Integer;
begin
  retval := xlBookSheetType(handle, index);
  result := SheetType(retval);
end;

function TXLBook.moveSheet(srcIndex: Integer; dstIndex: Integer): boolean;
begin
  result := xlBookMoveSheet(handle, srcIndex, dstIndex) > 0;
end;

function TXLBook.delSheet(index: Integer): boolean;
begin
  result := xlBookDelSheet(handle, index) > 0;
end;

function TXLBook.sheetCount: Integer;
begin
  result := xlBookSheetCount(handle);
end;

function TXLBook.addFormat: TXLFormat;
var
  fmtHandle: FormatHandle;
begin
  fmtHandle := xlBookAddFormat(handle, nil);
  if fmtHandle <> nil then
  begin
    result := TXLFormat.Create(fmtHandle, handle);
    formats := formats + [result];
  end
  else
    result := nil;
end;

function TXLBook.addFormat(initFormat: TXLFormat): TXLFormat;
var
  fmtHandle: FormatHandle;
begin
  fmtHandle := xlBookAddFormat(handle, initFormat.handle);
  if fmtHandle <> nil then
  begin
    result := TXLFormat.Create(fmtHandle, handle);
    formats := formats + [result];
  end
  else
    result := nil;
end;

function TXLBook.addFont: TXLFont;
var
  font: FontHandle;
begin
  font := xlBookAddFont(handle, nil);
  if font <> nil then
  begin
    result := TXLFont.Create(font, handle);
    fonts := fonts + [result];
  end
  else
    result := nil;
end;

function TXLBook.addFont(initFont: TXLFont): TXLFont;
var
  font: FontHandle;
begin
  font := xlBookAddFont(handle, initFont.handle);
  if font <> nil then
  begin
    result := TXLFont.Create(font, handle);
    fonts := fonts + [result];
  end
  else
    result := nil;
end;

function TXLBook.addRichString(): TXLRichString;
var
  richString: RichStringHandle;
begin
  richString := xlBookAddRichString(handle);
  if richString <> nil then
  begin
    result := TXLRichString.Create(richString, handle);
    richStrings := richStrings + [result];
  end
  else
    result := nil;
end;

function TXLBook.addCustomNumFormat(customNumFormat: PWideChar): Integer;
begin
  result := xlBookAddCustomNumFormat(handle, customNumFormat);
end;

function TXLBook.customNumFormat(fmt: Integer): WideString;
var  
  s: PWideChar;
begin
  s := xlBookCustomNumFormat(handle, fmt);
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

function TXLBook.format(index: Integer): TXLFormat;
var
  fmtHandle : FormatHandle;
begin
  fmtHandle := xlBookFormat(handle, index);
  if fmtHandle <> nil then
  begin
    result := TXLFormat.Create(fmtHandle, handle);
    formats := formats + [result];
  end
  else
    result := nil;
end;

function TXLBook.formatSize(): Integer;
begin
  result := xlBookFormatSize(handle);
end;

function TXLBook.font(index: Integer): TXLFont;
var
  fntHandle : FontHandle;
begin
  fntHandle := xlBookFont(handle, index);
  if fntHandle <> nil then
  begin
    result := TXLFont.Create(fntHandle, handle);
    fonts := fonts + [result];
  end
  else
   result := nil;
end;

function TXLBook.fontSize(): Integer;
begin
  result := xlBookFontSize(handle);
end;

function TXLBook.datePack(year, month, day: Integer): double;
begin
  result := xlBookDatePack(handle, year, month, day, 0, 0, 0, 0);
end;

function TXLBook.datePack(year, month, day, hour, min, sec: Integer): double;
begin
  result := xlBookDatePack(handle, year, month, day, hour, min, sec, 0);
end;

function TXLBook.datePack(year, month, day, hour, min, sec, msec: Integer): double;
begin
  result := xlBookDatePack(handle, year, month, day, hour, min, sec, msec);
end;

function TXLBook.dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer): boolean;
begin
  result := xlBookDateUnpack(handle, value, addr(year), addr(month), addr(day), nil, nil, nil, nil) > 0;
end;

function TXLBook.dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer; var hour: Integer; var min: Integer; var sec: Integer): boolean;
begin
  result := xlBookDateUnpack(handle, value, addr(year), addr(month), addr(day), addr(hour), addr(min), addr(sec), nil) > 0;
end;

function TXLBook.dateUnpack(value: double; var year: Integer; var month: Integer; var day: Integer; var hour: Integer; var min: Integer; var sec: Integer; var msec: Integer): boolean;
begin
  result := xlBookDateUnpack(handle, value, addr(year), addr(month), addr(day), addr(hour), addr(min), addr(sec), addr(msec)) > 0;
end;

function TXLBook.colorPack(red: Integer; green: Integer; blue: Integer): Color;
begin
  result := Color(xlBookColorPack(handle, red, green, blue));
end;

procedure TXLBook.colorUnpack(color: Color; var red: Integer; var green: Integer; var blue: Integer);
begin
  xlBookColorUnpack(handle, Integer(color), addr(red), addr(green), addr(blue));
end;

function TXLBook.activeSheet: Integer;
begin
  result := xlBookActiveSheet(handle);
end;

procedure TXLBook.setActiveSheet(index: Integer);
begin
  xlBookSetActiveSheet(handle, index);
end;

function TXLBook.pictureSize(): Integer;
begin
  result := xlBookPictureSize(handle);
end;

function TXLBook.getPicture(index: Integer; var Buffer: ByteArray; var size: Cardinal): PictureType;
var
  ptr: PByteArray;
  i: Integer;
  p: PByte;
  retval: Integer;
begin
  retval := xlBookGetPicture(handle, index, @ptr, @size);
  if retval <> $FF then
  begin
    SetLength(Buffer, size);
    p := @(ptr^);
    for i := 0 to size - 1 do
    begin
      Buffer[i] := Byte(p^);
      Inc(p);
    end;
  end;
  result := PictureType(retval);
end;

function TXLBook.addPicture(filename: PWideChar): Integer;
begin
  result := xlBookAddPicture(handle, filename);
end;

function TXLBook.addPicture2(data: PByteArray; size: Cardinal): Integer;
begin
  result := xlBookAddPicture2(handle, data, size);
end;

function TXLBook.addPictureAsLink(filename: PWideChar; insert: boolean): Integer;
begin
  result := xlBookAddPictureAsLink(handle, filename, Integer(insert));
end;

function TXLBook.defaultFont(var fontSize: Integer): WideString;
var
  s: PWideChar;
begin
  s := xlBookDefaultFont(handle, addr(fontSize));
  if s <> nil then
    result := WideString(s)
  else
    result := '';
end;

procedure TXLBook.setDefaultFont(fontName: PWideChar; fontSize: Integer);
begin
  xlBookSetDefaultFont(handle, fontName, fontSize);
end;

function TXLBook.GetRefR1C1: boolean;
begin
  result := (xlBookRefR1C1(handle) > 0);
end;

procedure TXLBook.SetRefR1C1(const value: boolean);
begin
  xlBookSetRefR1C1(handle, Integer(value));
end;

function TXLBook.GetDate1904: boolean;
begin
  result := (xlBookIsDate1904(handle) > 0);
end;

procedure TXLBook.SetDate1904(const value: boolean);
begin
  xlBookSetDate1904(handle, Integer(value));
end;

function TXLBook.GetTemplate: boolean;
begin
  result := (xlBookIsTemplate(handle) > 0);
end;

procedure TXLBook.SetTemplate(const value: boolean);
begin
  xlBookSetTemplate(handle, Integer(value));
end;

procedure TXLBook.setKey(name, key: PWideChar);
begin
  xlBookSetKey(handle, name, key);
end;

function TXLBook.rgbMode(): boolean;
begin
  result := boolean(xlBookRgbMode(handle));
end;

procedure TXLBook.setRgbMode(rgbMode: boolean);
begin
  xlBookSetRgbMode(handle, Integer(rgbMode));
end;

function TXLBook.calcMode(): CalcModeType;
begin
  result := CalcModeType(xlBookCalcMode(handle));
end;

procedure TXLBook.setCalcMode(calcMode: CalcModeType);
begin
  xlBookSetCalcMode(handle, Integer(calcMode));
end;

function TXLBook.biffVersion(): Integer;
begin
  result := xlBookBiffVersion(handle);
end;

function TXLBook.setLocale(locale: PAnsiChar): boolean;
begin
  result := xlBookSetLocale(handle, locale) > 0;
end;


function TXLBook.errorMessage: WideString;
begin
  result := WideString(xlBookErrorMessage(handle));
end;


procedure TXLBook.Release;
var
  i : Integer;
begin
  try
    for i:= 0 to High(fonts) do fonts[i].Destroy;
    for i:= 0 to High(formats) do formats[i].Destroy;
    for i:= 0 to High(richStrings) do richStrings[i].Destroy;
    for i:= 0 to High(sheets) do sheets[i].Destroy;
    xlBookRelease(handle);
  finally
  end;
end;

destructor TXLBook.Destroy;
begin
  Release;
  inherited Destroy;
end;

constructor TBinBook.Create();
begin
  inherited Create;
  handle := xlCreateBookC;
end;

constructor TXmlBook.Create();
begin
  inherited Create;
  handle := xlCreateXMLBookC;
end;


end.
