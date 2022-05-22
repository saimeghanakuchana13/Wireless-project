classdef NRUERx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;  % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grid b
        ofdmGrid;     % Before pre-coding nsc x nsym x nlayers
        
        % Channel and noise estimate
        noiseEst;
        chanEstGrid;
        % Channel estimation parameters
        sigFreq = 7;  % Channel smoothing in freq
        sigTime = 3;  % Channel smoothing in time
        lenFreq = 21;  % Filter length in freq
        Wtime; 
        
        % Recived symbols
        pdschChanEst;   % Channel estimate on the PDSCH
        pdschSymRaw;    % Raw symbols before equalization
        pdschSymEq;     % Equalized symbols
        dmrsSym;        % DM-RS Reference symbols
        
        % Slot number
        Nslot = 0;
        
        % RX beamforming vector.  
        rxBF;
        
        % Timing offset
        offset;
        
    end
    
    methods
        function obj = NRUERx(simParam, varargin)
            % Constructor
           
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the 
            % objects need to be converted to older structures.
            obj.carrierConfig = simParam.carrierConfig ;
            obj.pdschConfig = simParam.pdschConfig ;
            obj.waveformConfig = simParam.waveformConfig;
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end                        
            
        end

        function chanEst(obj, rxGrid)
            % Computes the channel estimate
            
            % TODO:  Get the TX DM-RS symbols and indices
            %   dmrsSymTx = ...
            %   dmrsInd = ...
            dmrsSymTx = nrPDSCHDMRS(obj.carrierConfig, obj.pdschConfig);
            dmrsInd = nrPDSCHDMRSIndices(obj.carrierConfig, obj.pdschConfig);
            
            % TODO:  Get RX symbols on the DM-RS
            %    dmrsSymRx = ...
            dmrsSymRx = rxGrid(dmrsInd);
            
            % TODO:  Get the raw channel estimate
            %   chanEstRaw = ...
            chanEstRaw = dmrsSymRx ./ dmrsSymTx;           
                        
            % Get the symbol numbers and sub-carrier indices of the
            % DM-RS symbols from the DM-RS
            %   dmrsSymNum(i) = symbol number for the i-th DM-RS symbol
            %   dmrsScInd(i) = sub-carrier index for the i-th DM-RS symbol
            [nsc, nsym, nlayers] = size(rxGrid);
            dmrsSymNum = floor(double(dmrsInd-1) / nsc)+1;
            dmrsScInd = mod((dmrsInd-1), nsc)+1;
            
            % TODO:  Get the list of all symbol numbers on which DM-RS was
            % transmitted.  You can use the unique command
            %   dmrsSymNums = unique(...);
            %   ndrmsSym = length(dmrsSymNums);
            dmrsSymNums = unique(dmrsSymNum);            
            ndrmsSym = length(dmrsSymNums);            

            % We first compute the channel and noise 
            % estimate on each of the symbols on which the DM-RS was 
            % transmitted.  We will store these in two arrays
            %   chanEstDmrs(k,i) = chan est on sub-carrier k in DM-RS
            %       symbol i
            %   noiseEstDmrs(i) = noise est for DM-RS symbol i
            chanEstDmrs = zeros(nsc, ndrmsSym);
            noiseEstDmrs  = zeros(ndrmsSym, 1);
            
            % Loop over the DM-RS symbols
            for i = 1:ndrmsSym
                
                % TODO:  Find the indices, k, in which the DM-RS
                % dmrsSymNum(k)= dmrsSymNum(i).
                %   I = find(...)                            
                I = find(dmrsSymNum == dmrsSymNums(i));  
                
                % TODO:  Get the sub-carrier indices and raw channel 
                % channel estimate for these RS on the symbol
                %   ind = ...
                %   raw = ...
                ind = dmrsScInd(I);
                raw = chanEstRaw(I);
                
                % TODO:  Use kernelReg to compute the channel estimate
                % on that DM-RS symbol.  Use the lenFreq and sigFreq
                % for the kernel length and sigma.
                %    chanEstDmrs(:,i) = kernelReg(...)
                chanEstDmrs(:,i) = kernelReg(ind, raw, nsc, ...
                    obj.lenFreq, obj.sigFreq);
                
                % TODO:  Compute the noise estimate on the symbol
                % using the residual method
                %    noiseEstDmrs(i) = ...
                noiseEstDmrs(i) = mean(abs(raw- chanEstDmrs(ind,i)).^2);
                
            end
            % TODO:  Find the noise estimate over the PDSCH by
            % averaging noiseEstDmrs
            %   obj.noiseEst = mean(...);         
            obj.noiseEst = mean(noiseEstDmrs);
                        
            % TODO:  Finally, we interpolate over time.
            % We will use an estimate of the form
            %    obj.chaneEstGrid = chanEstDrms*W
            % so that
            %    chanEstGrid(k,j) = \sum_i chanEstDmrs(k,i)*W(i,j)
            %
            % We use a kernel estimator
            %
            %     W(i,j) = W0(i,j) / \sum_k W0(k,j)
            %     W0(k,j) = exp(-D(k,j)^2/(2*obj.sigTime^2))
            %     D(k,j) = dmrsSymNum(k) - j
            %            
            D = dmrsSymNums - (1:nsym);
            W0 = exp(-D.^2/(2*obj.sigTime^2));
            W = W0 ./ sum(W0,1);
            
            % Save the time interpolation matrix
            obj.Wtime = W;                      
            
            % Create the channel estimate grid
            obj.chanEstGrid = chanEstDmrs*W;
            
        end
    end
    methods (Access = protected)
        
        function stepImpl(obj, y)
            
            % TODO:  Perform RX beamforming by multiplying y with the
            % the RX BF vector.  
            %   z = ...
            z = y*obj.rxBF;
                         
            % Get information for PDSCH and DM-RS allocations
            %[pdschIndices,dmrsIndices,dmrsSymbols,pdschIndicesInfo] = ...
            %    nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig); 
            [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(obj.carrierConfig, obj.pdschConfig);
            dmrsIndices = nrPDSCHDMRSIndices(obj.carrierConfig, obj.pdschConfig);
            dmrsSymbols = nrPDSCHDMRS(obj.carrierConfig, obj.pdschConfig);
            
            % Demodulate the RX signal
            obj.ofdmGrid = nrOFDMDemodulate(obj.carrierConfig, z);

            obj.chanEst(obj.ofdmGrid);
                  
            % Get channel estimate.
            % This is a poor channel estimate since we have not done
            % carrier and timing estimation.  But, this is OK for now.
            %[chanEstGrid, obj.noiseEst] = nrChannelEstimate(...
            %    obj.ofdmGrid,dmrsIndices,dmrsSymbols,...
             %   'CyclicPrefix',obj.carrierConfig.CyclicPrefix);  %,...
            
            % Extract raw symbols and channel estimate on PDSCH
            obj.pdschSymRaw = obj.ofdmGrid(pdschIndices);
            obj.pdschChanEst = obj.chanEstGrid(pdschIndices);
            obj.pdschSymEq = obj.pdschSymRaw ./ obj.pdschChanEst;
            
        end
        
    end
end

