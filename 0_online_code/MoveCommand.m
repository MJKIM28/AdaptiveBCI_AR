function MoveCommand(nx,ny)
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
% global mouse

mouse.mouseMove(0,0)
mouse.mouseMove(nx, ny); % BCI start button position