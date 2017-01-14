function  obj = Create_Object( T_list, ann_layers, bool_new_sample )
    
    obj.T_list = T_list;
    len = length(T_list);
    
    obj.net_list = cell(1, len);
    obj.pdf_list = cell(1,len);
    
    obj.ann_layers_list = ones(1,len).*ann_layers;
    
    obj.mean_list = zeros(1,len);
    obj.std_list = zeros(1,len);
    obj.hbw_list = zeros(1,len);
    
    obj.bool.new_sample = bool_new_sample;
    obj.bool.eval = 1;
    
end