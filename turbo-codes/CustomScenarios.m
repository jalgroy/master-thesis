% File CustomScenarios

% determine where your root directory is
load( 'CmlHome.mat' );

% determine where to store your files
base_name = 'Custom';
if ispc
    data_directory = strcat( '\output\', base_name, '\' );
else
    data_directory = strcat( '/output/', base_name, '/' );
end

full_directory = strcat( cml_home, data_directory );
if ~exist( full_directory, 'dir' )
    mkdir( full_directory);
end

record = 1;
sim_param(record).comment = 'Custom, r=1/2, K=128, S=4, iters=16';
sim_param(record).framesize = 128;
sim_param(record).sim_type = 'coded';
sim_param(record).code_configuration = 1;
sim_param(record).SNR = [-2:0.5:4];
sim_param(record).SNR_type = 'Eb/No in dB';
sim_param(record).modulation = 'BPSK';
sim_param(record).mod_order = 2;
sim_param(record).channel = 'AWGN';
sim_param(record).bicm = 1;
sim_param(record).demod_type = 0; 
sim_param(record).linetype = 'r-';
sim_param(record).legend = sim_param(record).comment;
sim_param(record).code_interleaver = ...
    strcat( 'CreateSRandomInterleaver(', int2str(sim_param(record).framesize ), ', 4)' );
sim_param(record).g1 = [1 0 0 1 1
    1 1 0 1 1];
sim_param(record).g2 = sim_param(record).g1;
sim_param(record).nsc_flag1 = 0;
sim_param(record).nsc_flag2 = 0;
sim_param(record).pun_pattern1 = [1 1
    1 0];
sim_param(record).pun_pattern2= [0 0
    0 1 ];
sim_param(record).tail_pattern1 = [1 1 1 1
    1 0 1 0];
sim_param(record).tail_pattern2 = [0 0 0 0
    0 1 0 1];
sim_param(record).max_iterations = 16;
sim_param(record).plot_iterations = sim_param(record).max_iterations;
sim_param(record).decoder_type = 0;
sim_param(record).filename = strcat( data_directory, 'CustomRate1by', ...
    int2str( floor( ( sum( sum( sim_param(record).pun_pattern1 ) ) + sum( sum( sim_param(record).pun_pattern2) ) )/size(sim_param(record).pun_pattern1,2) ) ), ...
    sim_param(record).channel, int2str( sim_param(record).framesize), '.mat' );
sim_param(record).reset = 0;;
sim_param(record).max_trials = 1e8*ones( size(sim_param(record).SNR) );
sim_param(record).minBER = 1e-5; 
sim_param(record).max_frame_errors = 10000*cat(2,[3300 100 4],ones( 1, length(sim_param(record).SNR)-3 ));
sim_param(record).save_rate = 50;

record = 2;
sim_param(record).comment = 'Custom, r=1/2, K=128, S=4, iters=16, PSTC v1';
sim_param(record).framesize = 128;
sim_param(record).sim_type = 'coded';
sim_param(record).code_configuration = 1;
sim_param(record).SNR = [-3:0.5:6];
sim_param(record).SNR_type = 'Eb/No in dB';
sim_param(record).modulation = 'BPSK';
sim_param(record).mod_order = 2;
sim_param(record).channel = 'AWGN';
sim_param(record).bicm = 1;
sim_param(record).demod_type = 0; 
sim_param(record).linetype = 'r-';
sim_param(record).legend = sim_param(record).comment;
sim_param(record).code_interleaver = ...
    strcat( 'CreateSRandomInterleaver(', int2str(sim_param(record).framesize ), ', 4)' );
sim_param(record).g1 = [1 0 1
    1 1 1];
sim_param(record).g2 = sim_param(record).g1;
sim_param(record).nsc_flag1 = 0;
sim_param(record).nsc_flag2 = 0;
sim_param(record).pun_pattern1 = [1 0 0 0 1 0 0 0
    1 1 1 1 0 1 1 1];
sim_param(record).pun_pattern2= [0 0 0 0 0 0 0 0
    0 1 1 1 1 1 1 1];
sim_param(record).tail_pattern1 = [1 1 1 1
    1 0 1 0];
sim_param(record).tail_pattern2 = [0 0 0 0
    0 1 0 1];
sim_param(record).max_iterations = 16;
sim_param(record).plot_iterations = sim_param(record).max_iterations;
sim_param(record).decoder_type = 0;
sim_param(record).filename = strcat( data_directory, 'Custom128Iter16NonSysV2Rate1by', ...
    int2str( floor( ( sum( sum( sim_param(record).pun_pattern1 ) ) + sum( sum( sim_param(record).pun_pattern2) ) )/size(sim_param(record).pun_pattern1,2) ) ), ...
    sim_param(record).channel, int2str( sim_param(record).framesize), '.mat' );
sim_param(record).reset = 0;;
sim_param(record).max_trials = 1e8*ones( size(sim_param(record).SNR) );
sim_param(record).minBER = 1e-5; 
sim_param(record).max_frame_errors = 1000*cat(2,[100 100 1000 1520 80 10],ones( 1, length(sim_param(record).SNR)-6 ));
sim_param(record).save_rate = 50;
