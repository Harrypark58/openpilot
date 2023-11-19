#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6650778850573399771);
void live_err_fun(double *nom_x, double *delta_x, double *out_938265339617516059);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_5640736375263316126);
void live_H_mod_fun(double *state, double *out_6262099730766365750);
void live_f_fun(double *state, double dt, double *out_4121089057125326820);
void live_F_fun(double *state, double dt, double *out_6765403946009136729);
void live_h_4(double *state, double *unused, double *out_7771952879335384996);
void live_H_4(double *state, double *unused, double *out_7268534510949444888);
void live_h_9(double *state, double *unused, double *out_4333434609376227545);
void live_H_9(double *state, double *unused, double *out_18684424315002582);
void live_h_10(double *state, double *unused, double *out_828795274476811386);
void live_H_10(double *state, double *unused, double *out_8061836490020891849);
void live_h_12(double *state, double *unused, double *out_8264527800838387238);
void live_H_12(double *state, double *unused, double *out_2249078102917483093);
void live_h_31(double *state, double *unused, double *out_4609540389192074352);
void live_H_31(double *state, double *unused, double *out_496484929407530616);
void live_h_32(double *state, double *unused, double *out_2922395548565996430);
void live_H_32(double *state, double *unused, double *out_5546711963287737114);
void live_h_13(double *state, double *unused, double *out_1750755361865855936);
void live_H_13(double *state, double *unused, double *out_2703051065181853727);
void live_h_14(double *state, double *unused, double *out_4333434609376227545);
void live_H_14(double *state, double *unused, double *out_18684424315002582);
void live_h_33(double *state, double *unused, double *out_8716799056665959868);
void live_H_33(double *state, double *unused, double *out_3647041934046388220);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}