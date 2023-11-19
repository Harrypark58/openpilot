#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_7139805516839873466);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4259303686840623864);
void car_H_mod_fun(double *state, double *out_8242501683440944479);
void car_f_fun(double *state, double dt, double *out_5097035619606452349);
void car_F_fun(double *state, double dt, double *out_3303120932540433895);
void car_h_25(double *state, double *unused, double *out_6726653390665798571);
void car_H_25(double *state, double *unused, double *out_3457977846036757342);
void car_h_24(double *state, double *unused, double *out_9118127512919358450);
void car_H_24(double *state, double *unused, double *out_1232270062057888780);
void car_h_30(double *state, double *unused, double *out_8818341386714792928);
void car_H_30(double *state, double *unused, double *out_3587316793179997412);
void car_h_26(double *state, double *unused, double *out_6993611477667413289);
void car_H_26(double *state, double *unused, double *out_7199481164910813566);
void car_h_27(double *state, double *unused, double *out_9115685362235055394);
void car_H_27(double *state, double *unused, double *out_5762080104980422323);
void car_h_29(double *state, double *unused, double *out_2136703353423548726);
void car_H_29(double *state, double *unused, double *out_3077085448865605228);
void car_h_28(double *state, double *unused, double *out_6495016988898430102);
void car_H_28(double *state, double *unused, double *out_5511812560284647105);
void car_h_31(double *state, double *unused, double *out_3838818032289703860);
void car_H_31(double *state, double *unused, double *out_3427331884159796914);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}