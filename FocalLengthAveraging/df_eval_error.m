function [df]=df_eval_error(f1,fgt)
%return error in focal length estimation,given estimated and ground truth
%values
err= f1/fgt-1;
df=abs(err);
end