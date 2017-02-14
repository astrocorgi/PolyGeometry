%This script takes the outputs from arolla_interp_poly.m and runs them
%through the addDebris function.

clear all 
close all
load polyA_interp

[interp_debrisa,interp_debrisb] = addDebris(polyA_interp,polyB_interp,1000,15,'interp_arolla.poly');

