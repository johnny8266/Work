#include "calculation.h"

calculation::calculation()
{
  //  set_add();
}


void calculation::set_add(double a)
{
  cal_a = 100. + a;
}


double calculation::get_add()
{
  return cal_a;
}
