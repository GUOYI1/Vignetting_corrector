
/* Header for color correction, created by Xie Donghai, 2015.12.14

*/



#ifndef COLOR_CORRECTION_HPP
#define COLOR_CORRECTION_HPP

#include <vector>
using namespace std;



_declspec(dllexport) int VignettingCorrectionUsingRG(unsigned char* pImage, int ht, int wd, vector<double>& vp);





#endif

