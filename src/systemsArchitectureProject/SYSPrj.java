package systemsArchitectureProject;

//imports for swing GUI
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

//imports for function and textfield class
//
//
//

public class SYSPrj implements ActionListener{
	
	// Creation of JFrame, TextField, Buttons and Button holders, and JPanel
	JFrame frame;
	JTextField textfield;
	JButton[] numberButtons = new JButton[10];
	JButton[] functionButtons = new JButton[3];
	JButton decButton, equButton, clrButton;
	JPanel panel;
	JPanel functionPanel;
	
	// Two fonts that are implemented in the Buttons and TextField
	Font myFont = new Font("Calbri", Font.BOLD,30);
	Font specialFont = new Font("Calbri", Font.BOLD,15);
	Font textFont = new Font("Calbri", Font.PLAIN, 25);
	
	// New Colors 
	Color LightRed = new Color(230, 116, 116);

	
	// Creation of class instance
	SYSPrj(){
		
		// Configuration for JFrame
		frame = new JFrame("SYSPrj");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(800, 680);
		frame.setLayout(null);
		frame.getContentPane().setBackground(Color.DARK_GRAY);
		
		// Configuration for back end TextField 
		textfield = new JTextField();
		textfield.setBounds(33, 25, 720, 350);
		textfield.setFont(textFont);
		
		// Special Button declaration for operations
		decButton = new JButton(".");
		equButton = new JButton("=");
		clrButton = new JButton("Clear");
		
		// Adding the function buttons to the array
		functionButtons[0] = decButton;
		functionButtons[1] = equButton;
		functionButtons[2] = clrButton;
		
		

		// Loop that formats the function buttons
		for(int i = 0; i < 3; i++) {
			functionButtons[i].addActionListener(this);
			functionButtons[i].setFont(specialFont);
			functionButtons[i].setFocusable(false);
			functionButtons[i].setBackground(Color.GRAY);
			functionButtons[i].setOpaque(true);
			functionButtons[i].setBorderPainted(false);
		}
		
		// Loop that formats the numerical buttons
		for(int i = 0; i < 10; i++) {
			numberButtons[i] = new JButton(String.valueOf(i));  //Sets Value of button to loop iteration
			numberButtons[i].addActionListener(this);
			numberButtons[i].setFont(myFont);
			numberButtons[i].setFocusable(false);
			numberButtons[i].setBackground(Color.LIGHT_GRAY);
		}
		
		// Creation of JPanel for the number buttons and operations
		panel = new JPanel();
		panel.setBounds(33,400,380,200);
		panel.setLayout(new GridLayout(3,4,8,8));
		panel.setBackground(Color.DARK_GRAY);
		
		// JPanel for functions and special operations (More will go here later
		functionPanel = new JPanel();
		functionPanel.setBounds(435, 400, 315, 200);
		functionPanel.setLayout(new GridLayout(4,4,20,20));
		functionPanel.setBackground(Color.LIGHT_GRAY);	
		
		
		// Placing all of the buttons on the JPanel
		panel.add(numberButtons[1]);
		panel.add(numberButtons[2]);
		panel.add(numberButtons[3]);
		panel.add(numberButtons[4]);
		panel.add(numberButtons[5]);
		panel.add(numberButtons[6]);
		panel.add(numberButtons[7]);
		panel.add(numberButtons[8]);
		panel.add(numberButtons[9]);
		panel.add(numberButtons[0]);
		panel.add(decButton);
		panel.add(equButton);
		
		// This is where function buttons will live
		functionPanel.add(clrButton);
		
		// Color formatting for buttons to make them stand out
		clrButton.setBackground(LightRed);
		equButton.setBackground(LightRed);
		
		// Adding panels to frame and setting frame to Visible
		frame.add(panel);
		frame.add(functionPanel);
		frame.add(textfield);
		frame.setVisible(true);		
	}
	
	// Main function
	public static void main(String[] args) {
			// Calling class 
			new SYSPrj();
	}
		
	
	// Listens for event to perform
	public void actionPerformed(ActionEvent e) {
		
		// Adds numbers to text field 
		for(int i = 0; i < 10; i++) {
			if(e.getSource() == numberButtons[i]) {
				textfield.setText(textfield.getText().concat(String.valueOf(i)));
			}
		}

		// Adds a decimal to the text field
		if(e.getSource() == decButton) {			
			textfield.setText(textfield.getText().concat("."));
		}

		// Clears text field if clear button is pressed
		if(e.getSource() == clrButton) {
			textfield.setText("");
		}
		
		// If equal button is pressed, goes through cases to determine operation
		if(e.getSource() == equButton) {	
			// Does Nothing rn
		}
	}
}
