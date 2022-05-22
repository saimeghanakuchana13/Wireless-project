classdef NRgNBTx < matlab.System
    % 5G NR gNB transmitter class    
    properties                
        carrierConfig;  % Carrier configuration         
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
  
        % Coded bits transmitted on PDSCH      
        txBits;         % Cell array of TX bits
        
        % OFDM grids
        ofdmGridLayer;     % Before pre-coding nsc x nsym x nlayers       
        ofdmGridAnt;       % Before pre-coding nsc x nsym x nantennas        
        
        % OFDM grid to visualize the type of symbols
        ofdmGridChan;
               
        % Transmitted data in last slots
        bits;           % TX bits
        pdschSym;       % TX symbols
        dmrsSym;        % TX data symbols
        
        % Channel
        chanNames;
        
        % Slot number
        Nslot = 0;
        
        % TX beamforming vector.  This is fixed.
        txBF;
        
        % HARQ Process 
        %nharq = 8;      % number of HARQ processes 
            
    end
    
    properties (Constant)
        % Indices for ofdmGridChan indicating the type of symbol
        
        
    end
    
    methods
        function obj = NRgNBTx(simParam, varargin)
            % Constructor
           
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the 
            % objects need to be converted to older structures.
            obj.carrierConfig = simParam.carrierConfig;
            obj.pdschConfig = simParam.pdschConfig;
            obj.waveformConfig =  simParam.waveformConfig;
                                    
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end   
      function setAck(obj, iharq)
            % Set that the HARQ transmission was received correctly
            obj.newDataAvail(iharq) = 1;                        
      end
    end
    methods (Access = protected)
    
        function x = stepImpl(obj)
            % step implementation.  Creates one slot of samples
            
            % Set the slot number if the PDSCH config
            obj.carrierConfig.NSlot = obj.Nslot;

            % Create the PDSCH grid before pre-coding
            nscPerRB = 12;
            nsc = obj.carrierConfig.NSizeGrid * nscPerRB;
            obj.ofdmGridLayer = zeros(nsc, ...
                obj.waveformConfig.SymbolsPerSlot, ...
                obj.pdschConfig.NumLayers);
            obj.ofdmGridChan = zeros(nsc, ...
                obj.waveformConfig.SymbolsPerSlot, ...
                obj.pdschConfig.NumLayers);
            
            % Get information for PDSCH and DM-RS allocations
            pdschIndices = nrPDSCHIndices(obj.carrierConfig, obj.pdschConfig);
            dmrsIndices = nrPDSCHDMRSIndices(obj.carrierConfig, obj.pdschConfig);
            obj.dmrsSym = nrPDSCHDMRS(obj.carrierConfig, obj.pdschConfig);
            
            % Get the PT-RS symbols and indices and insert them
            % in the TX grid
            ptrsSym = nrPDSCHPTRS(obj.carrierConfig, obj.pdschConfig);
            ptrsInd = nrPDSCHPTRSIndices(obj.carrierConfig, obj.pdschConfig);
            % Generate random bits
            bitsPerSym = 2;
            nsym = length(pdschIndices);
            nbits = bitsPerSym * nsym;
            obj.bits = randi([0 1], nbits, 1);
            obj.txBits = obj.bits;

            % Modulate the bits to symbols
            %M = nr.NRConst.modOrder(obj.pdschConfig.Modulation);
            M = 2^bitsPerSym;
            obj.pdschSym = qammod(obj.bits,M,'InputType','bit',...
                'UnitAveragePower',true);
            
            % Map symbols to OFDM grid
            obj.ofdmGridLayer(pdschIndices) = obj.pdschSym;
            obj.ofdmGridLayer(dmrsIndices) = obj.dmrsSym;
            % Map PT-RS symbols and indices and insert them
            % in the OFDM grid
           
            obj.ofdmGridLayer(ptrsInd) = ptrsSym;

          
            % Fill the channel with labels of the channels.
            % This is just for visualization
            obj.ofdmGridChan(pdschIndices) = 1;
            obj.ofdmGridChan(dmrsIndices) = 2;  
            obj.ofdmGridChan(ptrsInd) = 3;
            obj.chanNames = {'Other', 'PDSCH', 'DM-RS', 'PT-RS'}; 
            % Perform the OFDM modulation
            xlayer = nrOFDMModulate(obj.carrierConfig, obj.ofdmGridLayer);
      
            % TX beamforming:  At this point, 
            % xlayer will be an nsamp x 1 vector.  Use the TX beamforming
            % vector, obj.txBF, to map this to a  nsamp x nant matrix
            % where nant is the number of TX antennas.
            
            x = xlayer*obj.txBF';
            % Increment the slot number
            obj.Nslot = obj.Nslot + 1;
            
            
        end        
      
    end
end

