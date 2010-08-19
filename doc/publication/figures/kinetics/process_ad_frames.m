load('ts01_675_sample_image_frames.mat')

addpath('../../../../src/visualize_cell_features/');
all_images{675}{51} = [];
write_montage_image_set(all_images{675},'overall_set.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);

sampled_all = cell(0);
count = 1;
for i = round(2:(50-2)/35:50)
    sampled_all{count} = all_images{675}{i};
    count = count + 1;
end
write_montage_image_set(sampled_all,'sampled_set.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);

assembly_set = cell(0);
count = 1;
for i = round(3:(13-2)/9:13)
    assembly_set{count} = all_images{675}{i};
    count = count + 1;
end

write_montage_image_set(assembly_set,'assembly_set.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);

assembly_set = cell(0);
count = 1;
for i = round(3:(13-2)/4:13)
    assembly_set{count} = all_images{675}{i};
    count = count + 1;
end

write_montage_image_set(assembly_set,'assembly_set_4.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);

disassembly_set = cell(0);
count = 1;
start = 51-13;
last = 51;
for i = round(start:(last-start)/9:last)
    disassembly_set{count} = all_images{675}{i};
    count = count + 1;
end

write_montage_image_set(disassembly_set,'disassembly_set.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);

disassembly_set = cell(0);
count = 1;
start = 51-13;
last = 51;
for i = round(start:(last-start)/4:last)
    disassembly_set{count} = all_images{675}{i};
    count = count + 1;
end

write_montage_image_set(disassembly_set,'disassembly_set_4.png','bar_size',5,'pixel_size',0.215051,'bar_position',2);