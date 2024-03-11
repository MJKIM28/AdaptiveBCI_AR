function runCamera()
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;

mouse.mouseMove(4300, 730); %  button position

mouse.mousePress(InputEvent.BUTTON1_MASK); % left button
mouse.mouseRelease(InputEvent.BUTTON1_MASK); % left button
