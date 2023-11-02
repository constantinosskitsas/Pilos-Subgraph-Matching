object Form1: TForm1
  Left = 397
  Top = 200
  Width = 949
  Height = 167
  Caption = 'Form1'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -17
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 144
  TextHeight = 20
  object GenerateButton: TButton
    Left = 74
    Top = 37
    Width = 115
    Height = 38
    Caption = 'Generate'
    TabOrder = 0
    OnClick = GenerateButtonClick
  end
  object ExtractButton: TButton
    Left = 209
    Top = 37
    Width = 116
    Height = 38
    Caption = 'Extract'
    TabOrder = 1
    OnClick = ExtractButtonClick
  end
  object EditButton: TButton
    Left = 345
    Top = 37
    Width = 115
    Height = 38
    Caption = 'Edit'
    TabOrder = 2
    OnClick = EditButtonClick
  end
  object FormatButton: TButton
    Left = 480
    Top = 37
    Width = 115
    Height = 38
    Caption = 'Format'
    TabOrder = 3
    OnClick = FormatButtonClick
  end
  object CustomButton: TButton
    Left = 615
    Top = 37
    Width = 116
    Height = 38
    Caption = 'Custom'
    TabOrder = 4
    OnClick = CustomButtonClick
  end
  object InvoiceButton: TButton
    Left = 751
    Top = 37
    Width = 115
    Height = 38
    Caption = 'Invoice'
    TabOrder = 5
    OnClick = InvoiceButtonClick
  end
end
