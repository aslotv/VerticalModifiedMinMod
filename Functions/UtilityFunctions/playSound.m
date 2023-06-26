function playSound(varargin)
%PLAYSOUND is a function that loads audio MAT files and plays the loaded
%sound.
%   This function takes one input, the name of the soundfile to load, and
%   produces no output; it plays the sound in the soundfile. 
% 
% Input: nameOfSound - the name of the audio MAT file: 
%                      - playSound() or playSound use default name 'handel'
%                      - playSound('handel') or playSound handel
%                      - playSound('chirp') or playSound chirp
%                      - playSound('gong') or playSound gong
%                      - playSound('laughter') or playSound laughter
%                      - playSound('splat') or playSound splat
%                      - playSound('train') or playSound train
%
% Authors: Anita Stene Loetvedt and Dag L. Aksnes. 

switch nargin
    case 0
        nameOfSound = 'handel';
    case 1
        nameOfSound = varargin{1};
    otherwise
        error('Too many input arguments for playSound().')
end

load(nameOfSound,'y','Fs'); 
sound(y,Fs); 

end


