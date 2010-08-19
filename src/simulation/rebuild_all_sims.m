for i=10:20
    output_dir = fullfile('..','..','data','simulation',['phases_',num2str(i)],'Images','Paxillin'); 
    build_phase_data('output_dir', output_dir, 'phase_length', i)
end

for i=1:10
    output_dir = fullfile('..','..','data','simulation',['moving_',num2str(i)],'Images','Paxillin'); 
    build_moving_data('output_dir', output_dir, 'speed', i)
end

output_dir = fullfile('..','..','data','simulation','moving_0_5','Images','Paxillin');
build_moving_data('output_dir', output_dir, 'speed', 0.5)

build_stationary_data('output_dir', fullfile('..','..','data','simulation','stationary','Images','Paxillin'))