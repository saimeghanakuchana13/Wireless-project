classdef MIMOMPChan < matlab.System
    % MIMOMPChan:  MIMO multi-path fading channel    
    properties 
        fsamp;   % Sample rate in Hz
        
        % TX and RX antenna arrays
        txArr, rxArr;
                
        % Path properties
        aoaAz, aoaEl; % Angles of arrival in degrees
        aodAz, aodEl; % Angles of departure in degrees
        gain;  % path gains in dB
        dly;   % delays in seconds
        dop;   % doppler shift of each path in Hz
        
        % Fractional delay object
        fracDly;
        
        % Initial set of phases for the next step call
        phaseInit;
        
        % Previous gains and steering vectors
        utx, urx;
        gainTx, gainRx;
        
                                
    end
    
    methods 
        function obj = MIMOMPChan(varargin)
            % Constructor:  
            % The syntax allows you to call the constructor with syntax of
            % the form:
            %
            %     chan = SISOMPChan('Prop1', Val1, 'Prop2', val2, ...);
            if nargin >= 1
                obj.set(varargin{:});
            end
            
        end
        
    end
    methods (Access = protected)
        function setupImpl(obj)
              % setup:  This is called before the first step.
              
              % Create a dsp.VariableFractionalDelay object 
              obj.fracDly = dsp.VariableFractionalDelay(...
                'InterpolationMethod', 'Farrow','FilterLength',8,...
                'FarrowSmallDelayAction','Use off-centered kernel',...
                'MaximumDelay', 1024);                           
        end
        
        function resetImpl(obj)
            % reset:  Called on the first step after reset or release.
            
            % Reset the fracDly object
            obj.fracDly.reset();
            
            % Initialize phases, phaseInit, to a row vector of 
            % dimension equal to the number of paths with uniform values 
            % from 0 to 2pi
            npath = length(obj.gain);
            obj.phaseInit = 2*pi*rand(1,npath);
        end
        
        function releaseImpl(obj)
            % release:  Called after the release method
            
            % Release the fracDly object
            obj.fracDly.release();
        end
        
        function y = stepImpl(obj, x)
            % step:  Run samples through the channel
            % The input, x, should be nsamp x nanttx, 
            
                        
            % The delay in samples
            dlySamp = obj.dly;
            
            % Getting the TX steering vectors and element gains
            % along the angles of departure using the           
            [obj.utx, obj.gainTx] = obj.txArr.step(obj.aoaAz, obj.aoaEl);            
            [obj.urx, obj.gainRx] = obj.rxArr.step(obj.aodAz, obj.aodEl); 

            % Computing the total gain along each path in linear scale
            gainLin = db2mag(obj.gain + obj.gainTx + obj.gainRx);
                        
            % Initialize variables    
            nsamp  = size(x,1);
            nantrx = size(obj.urx,1);
            npath = length(obj.dly);
            y = zeros(nsamp,nantrx);
            
            % Getting the Doppler shift of each path from the TX and RX
            obj.dop = obj.txArr.doppler(obj.aoaAz, obj.aoaEl)+ obj.rxArr.doppler(obj.aodAz, obj.aodEl);
            % Using the Doppler shifts, compute the phase rotations 
            % on each path.  Specifically, if nsamp = length(x), create a
            % (nsamp+1) x npath matrix 
            %     phase(i,k) = phase rotation on sample i and path k
            nsamp = length(x);
            t = (0:1:nsamp)';
            t = t/obj.fsamp;
            oneMat = ones(nsamp+1,1);
            phase = oneMat*obj.phaseInit - t*obj.dop*2*pi;
            % Save the final phase, phase(nsamp+1,:)
            % as phaseInit for the next step.
            obj.phaseInit = phase(nsamp+1,:);
            % Loop over the paths
            for ipath = 1:npath

                % Computing the transmitted signal, x, along the 
                % TX spatial signature for path ipath. 
                %   z = nsamp x 1 vetor
                z = x*obj.utx(:,ipath);


                % Delaying the path by the dlySamp(ipath) using the fractional delay object
                zdly = obj.fracDly(z,dlySamp(ipath));

                % Multiplying by the gain 
                zdly = zdly*gainLin(ipath);

                % Phase rotation
                z1 = exp(1i*phase(1:nsamp,ipath)).*zdly;

                % Multiply by the RX spatial signature and add to y
                y = y + z1*obj.urx(:,ipath)';
            end
          
        end
    end
end