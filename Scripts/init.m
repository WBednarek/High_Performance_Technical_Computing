function init(obj)
% This method is called after the object is created.
% OBJ is the device object.
% End of function definition - DO NOT EDIT

% Extract the interface object.
interface = get(obj, 'Interface');

fclose(interface);

% Configure the buffer size to allow for waveform transfer.
set(interface, 'InputBufferSize', 12e6);
set(interface, 'OutputBufferSize', 12e6); % Originally is set to 50,000

expExplicitPoints10000 =[];