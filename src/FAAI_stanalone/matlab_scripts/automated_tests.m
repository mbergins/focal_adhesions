clc;
clear all;

addpath('matlab_scripts');

%Remove Edge Adhesions test

final_adhesion_goal = zeros(50);
final_adhesion_goal(20:30,20:30) = 1;

edge_adhesions_test = zeros(50);
edge_adhesions_test(1:10,1:10) = 1;
edge_adhesions_test(25:35,1:10) = 1;
edge_adhesions_test(40:50,1:10) = 1;

edge_adhesions_test(1:10,20:30) = 1;
edge_adhesions_test(20:30,20:30) = 1;
edge_adhesions_test(40:50,20:30) = 1;

edge_adhesions_test(1:10,40:50) = 1;
edge_adhesions_test(25:35,40:50) = 1;
edge_adhesions_test(40:50,40:50) = 1;

function_result = remove_edge_adhesions(edge_adhesions_test);

assert(all(all(final_adhesion_goal == function_result)),'Removing Edge Pixels: problem with removing edge pixels')

function_result = remove_edge_adhesions(final_adhesion_goal);
assert(all(all(final_adhesion_goal == function_result)),'Removing Edge Pixels: problem with removing non-edge adhesions')