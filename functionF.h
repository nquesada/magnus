/*! \file functionF.h
 * \brief Defines the  function that equals the product of the phase matching 
 * and the pump function. Note that here we have omitted the factor of 
 * of 1/\sqrt{pi} used in the main paper. That implies that the second order 
 * Magnus correction should be multiplied by 1/pi and the third order one by
 * 1/\sqrt{pi}^3. Note that since this function will be called many times
 * it should be implemented as efficiently as possible.
 * Input: the frequencies of the down converted modes, wo and we and the
 * frequency of the pump
 * Output: The product of the phase matching function and the pump function
 * at such frequencies.
 */


#ifndef _functionF_
#define _functionF_


double F(double wo, double we, double wp);
void c_F(double wo, double we, double wp, double res[2]);

#endif
