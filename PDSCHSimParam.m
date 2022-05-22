classdef PDSCHSimParam < matlab.mixin.SetGet
    % PDSCHSimParam:  Parameters for the PDSCHSimParam    
    properties
                
        % Settable parameters.  These should be set in the constructor
        NLayers = 1;    % number of layers
        Modulation = 'QPSK';  
        NRB = 51;  % number of resource blocks
        SubcarrierSpacing = 120;  % SCS in kHZ
        fc = 28e9;  % carrier frequency in Hz
        
        % Carrier parameters.  Right now, we assume only one BWP
        carrierConfig; 
                 
        % PDSCH parameters
        pdschConfig;
        
        % Waveform parameters
        waveformConfig;
         
    end
    
    methods
        function obj = PDSCHSimParam(varargin)
            % Constructor
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % Set the carrierParam settings
            obj.carrierConfig = nrCarrierConfig(...
                'NSizeGrid', obj.NRB, 'SubcarrierSpacing', obj.SubcarrierSpacing);
            %, ...'fc', obj.fc);
            
            % Compute the waveform parameters
            % We need to convert the carrier configuration to a structure
            obj.waveformConfig = nrOFDMInfo(obj.carrierConfig);
                                   
            % PDSCH parameters.  We allocate all the RBs and all 
            % the OFDM symbols in the slot
            res.Symbols = 0;  % Reserve first symbol
            res.PRB = (0:obj.NRB-1);      % Reserve all RBs
            res.Period = 1;  % Reserve on every slot
            obj.pdschConfig = nrPDSCHConfig(...
                'Modulation', 'QPSK', ...
                'PRBSet', (0:obj.NRB-1), ...
                'SymbolAllocation', [1,obj.waveformConfig.SymbolsPerSlot-1], ...
                'EnablePTRS', 1, 'PTRS', nrPDSCHPTRSConfig());
                %'Reserved', res);                                                     
        end
        
      
        
      
    end
end

