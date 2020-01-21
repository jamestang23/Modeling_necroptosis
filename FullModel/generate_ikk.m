function  v = generate_ikk(v)
if isfield(v,'IKK_input')

    time    = v.IKK_input{1};
    values  = v.IKK_input{2};    
           values1  = v.IKK_input{3};
else
    error('generate_ikk:generate_ikk() Incorrect IKK input trajs');
end
v.IKK_CURVE = interp1(time, values, 0:v.SIM_TIME,'pchip');
v.IKK_CURVE1 = interp1(time, values1, 0:v.SIM_TIME,'pchip');
end
