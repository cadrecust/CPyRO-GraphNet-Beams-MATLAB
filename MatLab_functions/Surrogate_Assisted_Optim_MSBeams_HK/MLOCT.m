function [Ao3model]=MLOCT(NLOCT,XTest,ATest,model2Use)
    
    % Model for prediction of As per cross-section
    if NLOCT==1
        if strcmp(model2Use,'CPyRO')
            % MLOCT 1
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_4000.mat");
            
            Ao3CPyRO1=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                        ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO1=extractdata(Ao3CPyRO1);
            
            Ao3model=[Ao3CPyRO1];
        else
            % MLOCT 1
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_4000.mat");
            paramGCNN=load("GCNN_As_Section_4000.mat");
            
            Ao3GNN1=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,...
                                            ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN1=extractdata(Ao3GNN1);
    
            Ao3model=[Ao3GNN1];
        end
    elseif NLOCT==2
        if strcmp(model2Use,'CPyRO')
            % MLOCT 1
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_MOConstrucT1_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_MOConstrucT1_4000.mat");
            
            Ao3CPyRO1=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                        ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO1=extractdata(Ao3CPyRO1);

            % MLOCT 2
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_MOConstrucT5_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_MOConstrucT5_4000.mat");
            
            Ao3CPyRO2=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                        ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO2=extractdata(Ao3CPyRO2);
            
            Ao3model=[Ao3CPyRO1,Ao3CPyRO2];
        else
            % MLOCT 1
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_MLOCT1_4000.mat");
            paramGCNN=load("GCNN_As_Section_MLOCT1_4000.mat");
            
            Ao3GNN1=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,...
                                            ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN1=extractdata(Ao3GNN1);
            
            % MLOCT 2
            
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_MLOCT5_4000.mat");
            paramGCNN=load("GCNN_As_Section_MLOCT5_4000.mat");
            
            Ao3GNN2=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,...
                                            ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN2=extractdata(Ao3GNN2);
            
            % Gather all estimations for each model
            Ao3model=[Ao3GNN1,Ao3GNN2];
        end
    elseif NLOCT==3
        if strcmp(model2Use,'CPyRO')
            % MLOCT 1
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_MOConstrucT1_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_MOConstrucT1_4000.mat");
            
            Ao3CPyRO1=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                    ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO1=extractdata(Ao3CPyRO1);

            % MLOCT 2
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_MOConstrucT3_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_MOConstrucT3_4000.mat");
            
            Ao3CPyRO2=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                        ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO2=extractdata(Ao3CPyRO2);

            % MLOCT 3
            nheadsparamnGATPIGNN=load("nHeads_GAT_PIGNN_As_Section_MOConstrucT5_4000.mat");
            paramPIGCNN=load("PIGCNN_As_Section_MOConstrucT5_4000.mat");
            
            Ao3CPyRO3=PIGNNmodel1fc1GAT1Conv1fc(paramPIGCNN.pignn,XTest,...
                                        ATest,nheadsparamnGATPIGNN.numHeads);
            Ao3CPyRO3=extractdata(Ao3CPyRO3);
            
            Ao3model=[Ao3CPyRO1,Ao3CPyRO2,Ao3CPyRO3];
        else
            % MLOCT 1
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_MLOCT1_4000.mat");
            paramGCNN=load("GCNN_As_Section_MLOCT1_4000.mat");
            
            Ao3GNN1=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,...
                                        ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN1=extractdata(Ao3GNN1);
            
            % MLOCT 2
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_MLOCT3_4000.mat");
            paramGCNN=load("GCNN_As_Section_MLOCT3_4000.mat");
            
            Ao3GNN2=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,XTest,...
                                            ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN2=extractdata(Ao3GNN2);

            % MLOCT 3
            nheadsparamnGATGNN=load("nHeads_GAT_GCNN_As_Section_MLOCT5_4000.mat");
            paramGCNN=load("GCNN_As_Section_MLOCT5_4000.mat");
            
            Ao3GNN3=PlainGNNmodel1fc1GAT1Conv1fc(paramGCNN.parameters,...
                                    XTest,ATest,nheadsparamnGATGNN.numHeads);
            Ao3GNN3=extractdata(Ao3GNN3);
            
            % Gather all estimations for each model
            Ao3model=[Ao3GNN1,Ao3GNN2,Ao3GNN3];
        end
    end
end