#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7139805516839873466) {
   out_7139805516839873466[0] = delta_x[0] + nom_x[0];
   out_7139805516839873466[1] = delta_x[1] + nom_x[1];
   out_7139805516839873466[2] = delta_x[2] + nom_x[2];
   out_7139805516839873466[3] = delta_x[3] + nom_x[3];
   out_7139805516839873466[4] = delta_x[4] + nom_x[4];
   out_7139805516839873466[5] = delta_x[5] + nom_x[5];
   out_7139805516839873466[6] = delta_x[6] + nom_x[6];
   out_7139805516839873466[7] = delta_x[7] + nom_x[7];
   out_7139805516839873466[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4259303686840623864) {
   out_4259303686840623864[0] = -nom_x[0] + true_x[0];
   out_4259303686840623864[1] = -nom_x[1] + true_x[1];
   out_4259303686840623864[2] = -nom_x[2] + true_x[2];
   out_4259303686840623864[3] = -nom_x[3] + true_x[3];
   out_4259303686840623864[4] = -nom_x[4] + true_x[4];
   out_4259303686840623864[5] = -nom_x[5] + true_x[5];
   out_4259303686840623864[6] = -nom_x[6] + true_x[6];
   out_4259303686840623864[7] = -nom_x[7] + true_x[7];
   out_4259303686840623864[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_8242501683440944479) {
   out_8242501683440944479[0] = 1.0;
   out_8242501683440944479[1] = 0;
   out_8242501683440944479[2] = 0;
   out_8242501683440944479[3] = 0;
   out_8242501683440944479[4] = 0;
   out_8242501683440944479[5] = 0;
   out_8242501683440944479[6] = 0;
   out_8242501683440944479[7] = 0;
   out_8242501683440944479[8] = 0;
   out_8242501683440944479[9] = 0;
   out_8242501683440944479[10] = 1.0;
   out_8242501683440944479[11] = 0;
   out_8242501683440944479[12] = 0;
   out_8242501683440944479[13] = 0;
   out_8242501683440944479[14] = 0;
   out_8242501683440944479[15] = 0;
   out_8242501683440944479[16] = 0;
   out_8242501683440944479[17] = 0;
   out_8242501683440944479[18] = 0;
   out_8242501683440944479[19] = 0;
   out_8242501683440944479[20] = 1.0;
   out_8242501683440944479[21] = 0;
   out_8242501683440944479[22] = 0;
   out_8242501683440944479[23] = 0;
   out_8242501683440944479[24] = 0;
   out_8242501683440944479[25] = 0;
   out_8242501683440944479[26] = 0;
   out_8242501683440944479[27] = 0;
   out_8242501683440944479[28] = 0;
   out_8242501683440944479[29] = 0;
   out_8242501683440944479[30] = 1.0;
   out_8242501683440944479[31] = 0;
   out_8242501683440944479[32] = 0;
   out_8242501683440944479[33] = 0;
   out_8242501683440944479[34] = 0;
   out_8242501683440944479[35] = 0;
   out_8242501683440944479[36] = 0;
   out_8242501683440944479[37] = 0;
   out_8242501683440944479[38] = 0;
   out_8242501683440944479[39] = 0;
   out_8242501683440944479[40] = 1.0;
   out_8242501683440944479[41] = 0;
   out_8242501683440944479[42] = 0;
   out_8242501683440944479[43] = 0;
   out_8242501683440944479[44] = 0;
   out_8242501683440944479[45] = 0;
   out_8242501683440944479[46] = 0;
   out_8242501683440944479[47] = 0;
   out_8242501683440944479[48] = 0;
   out_8242501683440944479[49] = 0;
   out_8242501683440944479[50] = 1.0;
   out_8242501683440944479[51] = 0;
   out_8242501683440944479[52] = 0;
   out_8242501683440944479[53] = 0;
   out_8242501683440944479[54] = 0;
   out_8242501683440944479[55] = 0;
   out_8242501683440944479[56] = 0;
   out_8242501683440944479[57] = 0;
   out_8242501683440944479[58] = 0;
   out_8242501683440944479[59] = 0;
   out_8242501683440944479[60] = 1.0;
   out_8242501683440944479[61] = 0;
   out_8242501683440944479[62] = 0;
   out_8242501683440944479[63] = 0;
   out_8242501683440944479[64] = 0;
   out_8242501683440944479[65] = 0;
   out_8242501683440944479[66] = 0;
   out_8242501683440944479[67] = 0;
   out_8242501683440944479[68] = 0;
   out_8242501683440944479[69] = 0;
   out_8242501683440944479[70] = 1.0;
   out_8242501683440944479[71] = 0;
   out_8242501683440944479[72] = 0;
   out_8242501683440944479[73] = 0;
   out_8242501683440944479[74] = 0;
   out_8242501683440944479[75] = 0;
   out_8242501683440944479[76] = 0;
   out_8242501683440944479[77] = 0;
   out_8242501683440944479[78] = 0;
   out_8242501683440944479[79] = 0;
   out_8242501683440944479[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5097035619606452349) {
   out_5097035619606452349[0] = state[0];
   out_5097035619606452349[1] = state[1];
   out_5097035619606452349[2] = state[2];
   out_5097035619606452349[3] = state[3];
   out_5097035619606452349[4] = state[4];
   out_5097035619606452349[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5097035619606452349[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5097035619606452349[7] = state[7];
   out_5097035619606452349[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3303120932540433895) {
   out_3303120932540433895[0] = 1;
   out_3303120932540433895[1] = 0;
   out_3303120932540433895[2] = 0;
   out_3303120932540433895[3] = 0;
   out_3303120932540433895[4] = 0;
   out_3303120932540433895[5] = 0;
   out_3303120932540433895[6] = 0;
   out_3303120932540433895[7] = 0;
   out_3303120932540433895[8] = 0;
   out_3303120932540433895[9] = 0;
   out_3303120932540433895[10] = 1;
   out_3303120932540433895[11] = 0;
   out_3303120932540433895[12] = 0;
   out_3303120932540433895[13] = 0;
   out_3303120932540433895[14] = 0;
   out_3303120932540433895[15] = 0;
   out_3303120932540433895[16] = 0;
   out_3303120932540433895[17] = 0;
   out_3303120932540433895[18] = 0;
   out_3303120932540433895[19] = 0;
   out_3303120932540433895[20] = 1;
   out_3303120932540433895[21] = 0;
   out_3303120932540433895[22] = 0;
   out_3303120932540433895[23] = 0;
   out_3303120932540433895[24] = 0;
   out_3303120932540433895[25] = 0;
   out_3303120932540433895[26] = 0;
   out_3303120932540433895[27] = 0;
   out_3303120932540433895[28] = 0;
   out_3303120932540433895[29] = 0;
   out_3303120932540433895[30] = 1;
   out_3303120932540433895[31] = 0;
   out_3303120932540433895[32] = 0;
   out_3303120932540433895[33] = 0;
   out_3303120932540433895[34] = 0;
   out_3303120932540433895[35] = 0;
   out_3303120932540433895[36] = 0;
   out_3303120932540433895[37] = 0;
   out_3303120932540433895[38] = 0;
   out_3303120932540433895[39] = 0;
   out_3303120932540433895[40] = 1;
   out_3303120932540433895[41] = 0;
   out_3303120932540433895[42] = 0;
   out_3303120932540433895[43] = 0;
   out_3303120932540433895[44] = 0;
   out_3303120932540433895[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3303120932540433895[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3303120932540433895[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3303120932540433895[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3303120932540433895[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3303120932540433895[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3303120932540433895[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3303120932540433895[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3303120932540433895[53] = -9.8000000000000007*dt;
   out_3303120932540433895[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3303120932540433895[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3303120932540433895[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3303120932540433895[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3303120932540433895[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3303120932540433895[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3303120932540433895[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3303120932540433895[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3303120932540433895[62] = 0;
   out_3303120932540433895[63] = 0;
   out_3303120932540433895[64] = 0;
   out_3303120932540433895[65] = 0;
   out_3303120932540433895[66] = 0;
   out_3303120932540433895[67] = 0;
   out_3303120932540433895[68] = 0;
   out_3303120932540433895[69] = 0;
   out_3303120932540433895[70] = 1;
   out_3303120932540433895[71] = 0;
   out_3303120932540433895[72] = 0;
   out_3303120932540433895[73] = 0;
   out_3303120932540433895[74] = 0;
   out_3303120932540433895[75] = 0;
   out_3303120932540433895[76] = 0;
   out_3303120932540433895[77] = 0;
   out_3303120932540433895[78] = 0;
   out_3303120932540433895[79] = 0;
   out_3303120932540433895[80] = 1;
}
void h_25(double *state, double *unused, double *out_6726653390665798571) {
   out_6726653390665798571[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3457977846036757342) {
   out_3457977846036757342[0] = 0;
   out_3457977846036757342[1] = 0;
   out_3457977846036757342[2] = 0;
   out_3457977846036757342[3] = 0;
   out_3457977846036757342[4] = 0;
   out_3457977846036757342[5] = 0;
   out_3457977846036757342[6] = 1;
   out_3457977846036757342[7] = 0;
   out_3457977846036757342[8] = 0;
}
void h_24(double *state, double *unused, double *out_9118127512919358450) {
   out_9118127512919358450[0] = state[4];
   out_9118127512919358450[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1232270062057888780) {
   out_1232270062057888780[0] = 0;
   out_1232270062057888780[1] = 0;
   out_1232270062057888780[2] = 0;
   out_1232270062057888780[3] = 0;
   out_1232270062057888780[4] = 1;
   out_1232270062057888780[5] = 0;
   out_1232270062057888780[6] = 0;
   out_1232270062057888780[7] = 0;
   out_1232270062057888780[8] = 0;
   out_1232270062057888780[9] = 0;
   out_1232270062057888780[10] = 0;
   out_1232270062057888780[11] = 0;
   out_1232270062057888780[12] = 0;
   out_1232270062057888780[13] = 0;
   out_1232270062057888780[14] = 1;
   out_1232270062057888780[15] = 0;
   out_1232270062057888780[16] = 0;
   out_1232270062057888780[17] = 0;
}
void h_30(double *state, double *unused, double *out_8818341386714792928) {
   out_8818341386714792928[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3587316793179997412) {
   out_3587316793179997412[0] = 0;
   out_3587316793179997412[1] = 0;
   out_3587316793179997412[2] = 0;
   out_3587316793179997412[3] = 0;
   out_3587316793179997412[4] = 1;
   out_3587316793179997412[5] = 0;
   out_3587316793179997412[6] = 0;
   out_3587316793179997412[7] = 0;
   out_3587316793179997412[8] = 0;
}
void h_26(double *state, double *unused, double *out_6993611477667413289) {
   out_6993611477667413289[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7199481164910813566) {
   out_7199481164910813566[0] = 0;
   out_7199481164910813566[1] = 0;
   out_7199481164910813566[2] = 0;
   out_7199481164910813566[3] = 0;
   out_7199481164910813566[4] = 0;
   out_7199481164910813566[5] = 0;
   out_7199481164910813566[6] = 0;
   out_7199481164910813566[7] = 1;
   out_7199481164910813566[8] = 0;
}
void h_27(double *state, double *unused, double *out_9115685362235055394) {
   out_9115685362235055394[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5762080104980422323) {
   out_5762080104980422323[0] = 0;
   out_5762080104980422323[1] = 0;
   out_5762080104980422323[2] = 0;
   out_5762080104980422323[3] = 1;
   out_5762080104980422323[4] = 0;
   out_5762080104980422323[5] = 0;
   out_5762080104980422323[6] = 0;
   out_5762080104980422323[7] = 0;
   out_5762080104980422323[8] = 0;
}
void h_29(double *state, double *unused, double *out_2136703353423548726) {
   out_2136703353423548726[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3077085448865605228) {
   out_3077085448865605228[0] = 0;
   out_3077085448865605228[1] = 1;
   out_3077085448865605228[2] = 0;
   out_3077085448865605228[3] = 0;
   out_3077085448865605228[4] = 0;
   out_3077085448865605228[5] = 0;
   out_3077085448865605228[6] = 0;
   out_3077085448865605228[7] = 0;
   out_3077085448865605228[8] = 0;
}
void h_28(double *state, double *unused, double *out_6495016988898430102) {
   out_6495016988898430102[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5511812560284647105) {
   out_5511812560284647105[0] = 1;
   out_5511812560284647105[1] = 0;
   out_5511812560284647105[2] = 0;
   out_5511812560284647105[3] = 0;
   out_5511812560284647105[4] = 0;
   out_5511812560284647105[5] = 0;
   out_5511812560284647105[6] = 0;
   out_5511812560284647105[7] = 0;
   out_5511812560284647105[8] = 0;
}
void h_31(double *state, double *unused, double *out_3838818032289703860) {
   out_3838818032289703860[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3427331884159796914) {
   out_3427331884159796914[0] = 0;
   out_3427331884159796914[1] = 0;
   out_3427331884159796914[2] = 0;
   out_3427331884159796914[3] = 0;
   out_3427331884159796914[4] = 0;
   out_3427331884159796914[5] = 0;
   out_3427331884159796914[6] = 0;
   out_3427331884159796914[7] = 0;
   out_3427331884159796914[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_7139805516839873466) {
  err_fun(nom_x, delta_x, out_7139805516839873466);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4259303686840623864) {
  inv_err_fun(nom_x, true_x, out_4259303686840623864);
}
void car_H_mod_fun(double *state, double *out_8242501683440944479) {
  H_mod_fun(state, out_8242501683440944479);
}
void car_f_fun(double *state, double dt, double *out_5097035619606452349) {
  f_fun(state,  dt, out_5097035619606452349);
}
void car_F_fun(double *state, double dt, double *out_3303120932540433895) {
  F_fun(state,  dt, out_3303120932540433895);
}
void car_h_25(double *state, double *unused, double *out_6726653390665798571) {
  h_25(state, unused, out_6726653390665798571);
}
void car_H_25(double *state, double *unused, double *out_3457977846036757342) {
  H_25(state, unused, out_3457977846036757342);
}
void car_h_24(double *state, double *unused, double *out_9118127512919358450) {
  h_24(state, unused, out_9118127512919358450);
}
void car_H_24(double *state, double *unused, double *out_1232270062057888780) {
  H_24(state, unused, out_1232270062057888780);
}
void car_h_30(double *state, double *unused, double *out_8818341386714792928) {
  h_30(state, unused, out_8818341386714792928);
}
void car_H_30(double *state, double *unused, double *out_3587316793179997412) {
  H_30(state, unused, out_3587316793179997412);
}
void car_h_26(double *state, double *unused, double *out_6993611477667413289) {
  h_26(state, unused, out_6993611477667413289);
}
void car_H_26(double *state, double *unused, double *out_7199481164910813566) {
  H_26(state, unused, out_7199481164910813566);
}
void car_h_27(double *state, double *unused, double *out_9115685362235055394) {
  h_27(state, unused, out_9115685362235055394);
}
void car_H_27(double *state, double *unused, double *out_5762080104980422323) {
  H_27(state, unused, out_5762080104980422323);
}
void car_h_29(double *state, double *unused, double *out_2136703353423548726) {
  h_29(state, unused, out_2136703353423548726);
}
void car_H_29(double *state, double *unused, double *out_3077085448865605228) {
  H_29(state, unused, out_3077085448865605228);
}
void car_h_28(double *state, double *unused, double *out_6495016988898430102) {
  h_28(state, unused, out_6495016988898430102);
}
void car_H_28(double *state, double *unused, double *out_5511812560284647105) {
  H_28(state, unused, out_5511812560284647105);
}
void car_h_31(double *state, double *unused, double *out_3838818032289703860) {
  h_31(state, unused, out_3838818032289703860);
}
void car_H_31(double *state, double *unused, double *out_3427331884159796914) {
  H_31(state, unused, out_3427331884159796914);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
