function [bin] = alignIntervalB1(interval,mzTest,SpTest, ppm, binOfLengthOneOK,splitCorrection) 

% Eppm is the ppm error window
Eppm = @(mz, ppm) (2*ppm*mz)/1000000;
% Initialise the bin that holds the alignment results
field1 = 'mz';  value1 = zeros(1,10);
field2 = 'spectra';  value2 =zeros(1,10);
field3 = 'centroidspectra';  value3 =zeros(1,1);
field4 = 'centroidmz';  value4 =zeros(1,1);
bin = struct(field1,value1,field2,value2, field3,value3,field4,value4)   ;  
counter_bin = 1;

% For all ions in the interval do the alignment
        if (size(interval.ions,1)>1)
            % Scan the interval(a part of reference spectrum)
            indx_temp = 1;
            ind_strt = max(indx_temp);
            indx_temp = find(interval.ions <= interval.ions(ind_strt)...
                + Eppm(interval.ions(ind_strt),ppm));
            % Ba is the first window 
            Ba = interval.ions(indx_temp);
            BaSpectra = interval.spectra(indx_temp);
            number_of_ions_in_intervali = length(interval.ions);
            
            while(ind_strt <= number_of_ions_in_intervali)
                        ind_strt = max(indx_temp)+1;
                        if(ind_strt <= number_of_ions_in_intervali)
                                            indx_temp = find(interval.ions(ind_strt:end)...
                                                <= interval.ions(ind_strt) + Eppm(interval.ions(ind_strt),ppm)) + (ind_strt-1);
                                            % Bb is the second window
                                            Bb = interval.ions(indx_temp);
                                            BbSpectra = interval.spectra(indx_temp);
                                            %check if the biggest ion in Ba
                                            %is closer than Eppm to
                                            %smallest ion in Bb if yes
                                            %overlap_flag=1
                                            overlap_flag = Bb(1) <= (Ba(end)+Eppm(Ba(end),ppm));
                                    %         4.Combine Ba and Bb into a bin B and perform
                                    %         hierarchical clustering on the peaks. One may
                                    %         get up to three clusters in this process returned
                                    %         in an array, C.
                                            if(overlap_flag)
                                                        %4.i clustering
                                                        Bab = [Ba ; Bb];
                                                        BabSpectra = [BaSpectra ; BbSpectra];
                                                        [~,D] = linkage4MStest(Bab,ppm);
                                                        numberOfClusters = numel(unique(D));
                                                        %4.ii. Bb ‹- Last cluster in C
                                                        Bb = Bab(D == D(end));
                                                        BbSpectra = BabSpectra(D == D(end));
                                                        %4.iii. If |C| == 2
                                                        %Ba ‹- C[2]
                                                        %4.v. If |C| == 3 print C[1] to output and Ba ‹- C[2]
                                                        switch numberOfClusters
                                                            case 1
                                                                Ba = [];
                                                                BaSpectra = [];
                                                            case 2
                                                                Ba = Bab(D == D(1));
                                                                BaSpectra = BabSpectra(D == D(1)); 
                                                            case 3
                                                                % writes
                                                                % the first
                                                                % cluster
                                                                % to the
                                                                % output
                                                               Ba = Bab(D == D(1));
                                                               BaSpectra = BabSpectra(D == D(1)); 
                                                               [temp1 temp2 temp3] = ind2sub(size(mzTest) , BaSpectra);
                                                               for (l = 1:length(Ba))
                                                                        temp_spectra(l) = SpTest(temp1(l), temp2(l), temp3(l));
                                                               end
                                                                centroidBasp = centroid(temp_spectra, Ba);
                                                                centroidBamz = centroid(Ba, temp_spectra);
                                                               % write to the bin                 
                                                               if(numel(Ba)>1)
                                                                    bin(counter_bin).mz = Ba;
                                                                    bin(counter_bin).spectra = BaSpectra;
                                                                    bin(counter_bin).centroidmz = centroidBamz ;
                                                                    bin(counter_bin).centroidspectra = centroidBasp ;
                                                                    counter_bin = counter_bin + 1;
                                                                    clear temp_spectra temp1 temp2 temp3 centroidBasp centroidBamz
                                                               else
                                                                   if(binOfLengthOneOK)
                                                                    bin(counter_bin).mz = Ba;
                                                                    bin(counter_bin).spectra = BaSpectra;
                                                                    bin(counter_bin).centroidmz = centroidBamz ;
                                                                    bin(counter_bin).centroidspectra = centroidBasp ;
                                                                    counter_bin = counter_bin + 1;
                                                                    clear temp_spectra temp1 temp2 temp3 centroidBasp centroidBamz

                                                                   end
                                                               end
                                                                       
                                                                % Ba = second last cluster in C    
                                                                Ba = Bab(find((D~= D(1) & D~=D(end))));
                                                                BaSpectra = BabSpectra(find((D~= D(1) & D~=D(end))));
                                                                
                                                                % recunstruct
                                                                % Bb
                                                                            
                                                                a_ind_strt = max(indx_temp)+1;
                                                                if(a_ind_strt <= number_of_ions_in_intervali  & interval.ions(a_ind_strt) < (Bb(1)+Eppm(Bb(1),ppm)))
                                                                         ind_strt = max(indx_temp)+1;
                                                                         indx_temp = find(interval.ions(ind_strt:end)...
                                                                        <= (Bb(1)+Eppm(Bb(1),ppm))) + (ind_strt-1);
                                                                    % Bb is the second window
                                                                    Bb = [Bb ; interval.ions(indx_temp)];
                                                                    BbSpectra = [BbSpectra; interval.spectra(indx_temp)];
                                                                end
                                                        end
                                                                                                 
                                            end
                                                        
                                           % write Ba to output(bin) 
                                           if(~isempty(Ba))
                                                        [temp1 temp2 temp3] = ind2sub(size(mzTest) , BaSpectra);
                                                        for (l = 1:length(Ba))
                                                                temp_spectra(l) = SpTest(temp1(l), temp2(l), temp3(l));
                                                        end
                                                        centroidBasp = centroid(temp_spectra, Ba);
                                                        centroidBamz = centroid(Ba, temp_spectra);
                                                        clear temp_spectra
                                                         

                                                               if(numel(Ba)>1)
                                                                    bin(counter_bin).mz = Ba;
                                                                    bin(counter_bin).spectra = BaSpectra;
                                                                    bin(counter_bin).centroidmz = centroidBamz ;
                                                                    bin(counter_bin).centroidspectra = centroidBasp ;
                                                                    counter_bin = counter_bin + 1;
                                                                    clear temp_spectra temp1 temp2 temp3 centroidBasp centroidBamz
                                                               else
                                                                   if(binOfLengthOneOK)
                                                                    bin(counter_bin).mz = Ba;
                                                                    bin(counter_bin).spectra = BaSpectra;
                                                                    bin(counter_bin).centroidmz = centroidBamz ;
                                                                    bin(counter_bin).centroidspectra = centroidBasp ;
                                                                    counter_bin = counter_bin + 1;
                                                                    clear temp_spectra temp1 temp2 temp3 centroidBasp centroidBamz

                                                                   end
                                                               end
                                           end
                                            % Move Ba one window up (Ba=Bb)
                                            Ba = Bb;
                                            BaSpectra = BbSpectra;
                        end
            end
                       
        else
            if(binOfLengthOneOK)
                 bin(counter_bin).mz = interval.ions;
                 bin(counter_bin).spectra = interval.spectra;
                 bin(counter_bin).centroidmz = interval.ions ;
                 bin(counter_bin).centroidspectra = SpTest(interval.I, interval.J, interval.Z) ;
                 counter_bin = counter_bin + 1;
            end
          end
                        
           
        
        if(splitCorrection)
        %% merge split bins
        % get bins centroids
        cent = cell2mat({bin.centroidmz});
        % find empty bins 
        ind_zero = find(cent==0);
        bin(ind_zero) = [];
        cent(ind_zero) = [];
        % find diftance between bins
        wdiff = diff(cent);
        D = wdiff;
        % scale the distance to ppm 
        for j=1:length(wdiff)
            D(j)= D(j)/Eppm(cent(j),ppm/2);
        end
        % find split bins
        ind = find(D<1);
        % merge bins
        for i=1:(numel(ind))
                bin(ind(i)).mz = [bin(ind(i)).mz; bin(ind(i)+1).mz];
                bin(ind(i)).spectra = [bin(ind(i)).spectra; bin(ind(i)+1).spectra];
                bin(ind(i)).centroidmz = centroid([bin(ind(i)).centroidmz bin(ind(i)+1).centroidmz],[bin(ind(i)).centroidspectra bin(ind(i)+1).centroidspectra]);
                bin(ind(i)).centroidspectra = centroid([bin(ind(i)).centroidspectra bin(ind(i)+1).centroidspectra],[bin(ind(i)).centroidmz bin(ind(i)+1).centroidmz]);

        end
        % delete the split part of the bin
        ind  = ind+1;
        bin(ind) = [];
        end
        %%
        clearvars -except bin
end