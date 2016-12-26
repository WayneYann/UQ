function net = train_ann_idt(X0, IDT0)

    %     X0 = dlmread('data/samples.txt');
    %     IDT0 = dlmread( 'data/samples_out_ign.txt' );
    %     IDT0 = log(IDT0);
    
    %Train ANN
    hiddenLayerSize = 13;
    net = feedforwardnet(hiddenLayerSize,'trainlm');
    net.trainParam.epochs = 500;
    net.trainParam.goal = 1e-6;
    net = train(net, X0', IDT0');
end