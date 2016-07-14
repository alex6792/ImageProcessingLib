/*!
 * \author Alexandre Krebs
 * \file color.hpp
 * \brief color class
 */


#pragma once

#ifndef COLOR_HPP
#define COLOR_HPP


#include <iostream>


/*!
 * \class Color
 * \brief class representing a color by the three channels RGB and the opacity channel alpha
 */
class Color
{
    private :
        unsigned char R; /*!< red component*/
        unsigned char G; /*!< green component*/
        unsigned char B; /*!< blue component*/
        unsigned char A; /*!< alpha component*/

    public :
        // constructors & accessors
        Color();
        Color(unsigned char, unsigned char, unsigned char);
        Color(unsigned char, unsigned char, unsigned char, unsigned char);
        unsigned char red() const;
        unsigned char green() const;
        unsigned char blue() const;
        unsigned char alpha() const;
        unsigned char gray() const;
        int hue() const;
        int saturation() const;
        unsigned char value() const;

        void red(unsigned char);
        void green(unsigned char);
        void blue(unsigned char);
        void alpha(unsigned char);
        void HSV2RGB(int, int, unsigned char);

        // operators
        void operator+=(const Color&);
        void operator-=(const Color&);
        void operator*=(const Color&);

        void operator+=(const unsigned char&);
        void operator-=(const unsigned char&);
        void operator*=(const unsigned char&);
        void operator/=(const unsigned char&);

        Color operator+(const Color&) const;
        Color operator-(const Color&) const;
        Color operator*(const Color&) const;

        Color operator+(const unsigned char&) const;
        Color operator-(const unsigned char&) const;
        Color operator*(const unsigned char&) const;
        Color operator/(const unsigned char&) const;

        bool operator==(const Color&) const;
        bool operator!=(const Color&) const;
};

// operators
Color operator+(const unsigned char&, const Color&);
Color operator-(const unsigned char&, const Color&);
Color operator*(const unsigned char&, const Color&);
std::ostream& operator<<(std::ostream&, const Color&);

// predefined colors
#define AliceBlue Color(0xF0,0xF8,0xFF)
#define AntiqueWhite Color(0xFA,0xEB,0xD7)
#define Aqua Color(0x00,0xFF,0xFF)
#define Aquamarine Color(0x7F,0xFF,0xD4)
#define Azure Color(0xF0,0xFF,0xFF)
#define Beige Color(0xF5,0xF5,0xDC)
#define Bisque Color(0xFF,0xE4,0xC4)
#define Black Color(0x00,0x00,0x00)
#define BlanchedAlmond Color(0xFF,0xEB,0xCD)
#define Blue Color(0x00,0x00,0xFF)
#define BlueViolet Color(0x8A,0x2B,0xE2)
#define Brown Color(0xA5,0x2A,0x2A)
#define BurlyWood Color(0xDE,0xB8,0x87)
#define CadetBlue Color(0x5F,0x9E,0xA0)
#define Chartreuse Color(0x7F,0xFF,0x00)
#define Chocolate Color(0xD2,0x69,0x1E)
#define Coral Color(0xFF,0x7F,0x50)
#define CornflowerBlue Color(0x64,0x95,0xED)
#define Cornsilk Color(0xFF,0xF8,0xDC)
#define Crimson Color(0xDC,0x14,0x3C)
#define Cyan Color(0x00,0xFF,0xFF)
#define DarkBlue Color(0x00,0x00,0x8B)
#define DarkCyan Color(0x00,0x8B,0x8B)
#define DarkGoldenRod Color(0xB8,0x86,0x0B)
#define DarkGray Color(0xA9,0xA9,0xA9)
#define DarkGreen Color(0x00,0x64,0x00)
#define DarkKhaki Color(0xBD,0xB7,0x6B)
#define DarkMagenta Color(0x8B,0x00,0x8B)
#define DarkOliveGreen Color(0x55,0x6B,0x2F)
#define DarkOrange Color(0xFF,0x8C,0x00)
#define DarkOrchid Color(0x99,0x32,0xCC)
#define DarkRed Color(0x8B,0x00,0x00)
#define DarkSalmon Color(0xE9,0x96,0x7A)
#define DarkSeaGreen Color(0x8F,0xBC,0x8F)
#define DarkSlateBlue Color(0x48,0x3D,0x8B)
#define DarkSlateGray Color(0x2F,0x4F,0x4F)
#define DarkTurquoise Color(0x00,0xCE,0xD1)
#define DarkViolet Color(0x94,0x00,0xD3)
#define DeepPink Color(0xFF,0x14,0x93)
#define DeepSkyBlue Color(0x00,0xBF,0xFF)
#define DimGray Color(0x69,0x69,0x69)
#define DodgerBlue Color(0x1E,0x90,0xFF)
#define FireBrick Color(0xB2,0x22,0x22)
#define FloralWhite Color(0xFF,0xFA,0xF0)
#define ForestGreen Color(0x22,0x8B,0x22)
#define Fuchsia Color(0xFF,0x00,0xFF)
#define Gainsboro Color(0xDC,0xDC,0xDC)
#define GhostWhite Color(0xF8,0xF8,0xFF)
#define Gold Color(0xFF,0xD7,0x00)
#define GoldenRod Color(0xDA,0xA5,0x20)
#define Gray Color(0x80,0x80,0x80)
#define Green Color(0x00,0x80,0x00)
#define GreenYellow Color(0xAD,0xFF,0x2F)
#define HoneyDew Color(0xF0,0xFF,0xF0)
#define HotPink Color(0xFF,0x69,0xB4)
#define IndianRed Color(0xCD,0x5C,0x5C)
#define Indigo Color(0x4B,0x00,0x82)
#define Ivory Color(0xFF,0xFF,0xF0)
#define Khaki Color(0xF0,0xE6,0x8C)
#define Lavender Color(0xE6,0xE6,0xFA)
#define LavenderBlush Color(0xFF,0xF0,0xF5)
#define LawnGreen Color(0x7C,0xFC,0x00)
#define LemonChiffon Color(0xFF,0xFA,0xCD)
#define LightBlue Color(0xAD,0xD8,0xE6)
#define LightCoral Color(0xF0,0x80,0x80)
#define LightCyan Color(0xE0,0xFF,0xFF)
#define LightGoldenRodYellow Color(0xFA,0xFA,0xD2)
#define LightGray Color(0xD3,0xD3,0xD3)
#define LightGreen Color(0x90,0xEE,0x90)
#define LightPink Color(0xFF,0xB6,0xC1)
#define LightSalmon Color(0xFF,0xA0,0x7A)
#define LightSeaGreen Color(0x20,0xB2,0xAA)
#define LightSkyBlue Color(0x87,0xCE,0xFA)
#define LightSlateGray Color(0x77,0x88,0x99)
#define LightSteelBlue Color(0xB0,0xC4,0xDE)
#define LightYellow Color(0xFF,0xFF,0xE0)
#define Lime Color(0x00,0xFF,0x00)
#define LimeGreen Color(0x32,0xCD,0x32)
#define Linen Color(0xFA,0xF0,0xE6)
#define Magenta Color(0xFF,0x00,0xFF)
#define Maroon Color(0x80,0x00,0x00)
#define MediumAquaMarine Color(0x66,0xCD,0xAA)
#define MediumBlue Color(0x00,0x00,0xCD)
#define MediumOrchid Color(0xBA,0x55,0xD3)
#define MediumPurple Color(0x93,0x70,0xDB)
#define MediumSeaGreen Color(0x3C,0xB3,0x71)
#define MediumSlateBlue Color(0x7B,0x68,0xEE)
#define MediumSpringGreen Color(0x00,0xFA,0x9A)
#define MediumTurquoise Color(0x48,0xD1,0xCC)
#define MediumVioletRed Color(0xC7,0x15,0x85)
#define MidnightBlue Color(0x19,0x19,0x70)
#define MintCream Color(0xF5,0xFF,0xFA)
#define MistyRose Color(0xFF,0xE4,0xE1)
#define Moccasin Color(0xFF,0xE4,0xB5)
#define NavajoWhite Color(0xFF,0xDE,0xAD)
#define Navy Color(0x00,0x00,0x80)
#define OldLace Color(0xFD,0xF5,0xE6)
#define Olive Color(0x80,0x80,0x00)
#define OliveDrab Color(0x6B,0x8E,0x23)
#define Orange Color(0xFF,0xA5,0x00)
#define OrangeRed Color(0xFF,0x45,0x00)
#define Orchid Color(0xDA,0x70,0xD6)
#define PaleGoldenRod Color(0xEE,0xE8,0xAA)
#define PaleGreen Color(0x98,0xFB,0x98)
#define PaleTurquoise Color(0xAF,0xEE,0xEE)
#define PaleVioletRed Color(0xDB,0x70,0x93)
#define PapayaWhip Color(0xFF,0xEF,0xD5)
#define PeachPuff Color(0xFF,0xDA,0xB9)
#define Peru Color(0xCD,0x85,0x3F)
#define Pink Color(0xFF,0xC0,0xCB)
#define Plum Color(0xDD,0xA0,0xDD)
#define PowderBlue Color(0xB0,0xE0,0xE6)
#define Purple Color(0x80,0x00,0x80)
#define RebeccaPurple Color(0x66,0x33,0x99)
#define Red Color(0xFF,0x00,0x00)
#define RosyBrown Color(0xBC,0x8F,0x8F)
#define RoyalBlue Color(0x41,0x69,0xE1)
#define SaddleBrown Color(0x8B,0x45,0x13)
#define Salmon Color(0xFA,0x80,0x72)
#define SandyBrown Color(0xF4,0xA4,0x60)
#define SeaGreen Color(0x2E,0x8B,0x57)
#define SeaShell Color(0xFF,0xF5,0xEE)
#define Sienna Color(0xA0,0x52,0x2D)
#define Silver Color(0xC0,0xC0,0xC0)
#define SkyBlue Color(0x87,0xCE,0xEB)
#define SlateBlue Color(0x6A,0x5A,0xCD)
#define SlateGray Color(0x70,0x80,0x90)
#define Snow Color(0xFF,0xFA,0xFA)
#define SpringGreen Color(0x00,0xFF,0x7F)
#define SteelBlue Color(0x46,0x82,0xB4)
#define Tan Color(0xD2,0xB4,0x8C)
#define Teal Color(0x00,0x80,0x80)
#define Thistle Color(0xD8,0xBF,0xD8)
#define Tomato Color(0xFF,0x63,0x47)
#define Turquoise Color(0x40,0xE0,0xD0)
#define Violet Color(0xEE,0x82,0xEE)
#define Wheat Color(0xF5,0xDE,0xB3)
#define White Color(0xFF,0xFF,0xFF)
#define WhiteSmoke Color(0xF5,0xF5,0xF5)
#define Yellow Color(0xFF,0xFF,0x00)
#define YellowGreen Color(0x9A,0xCD,0x32)


#endif // COLOR_HPP