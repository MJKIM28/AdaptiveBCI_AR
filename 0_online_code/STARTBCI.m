function STARTBCI()
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
% global mouse

mouse.mouseMove(0,0)
mouse.mouseMove(350, 780); % BCI start button position
mouse.mousePress(InputEvent.BUTTON1_MASK); % right button
mouse.mouseRelease(InputEvent.BUTTON1_MASK); % right button