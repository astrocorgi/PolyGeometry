%This script takes the outputs from arolla_interp_poly.m and runs them
%through the addDebris function.

load polyA_interp

addDebris(polyA_interp,polyB_interp,1000,15);