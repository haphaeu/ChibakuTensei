/*
 * Orbits
 *
 * Draws the orbits of planets around the sun.
 *
 * Calculations are done using classic Newton mechanics with planet data
 * from NASA
 * https://nssdc.gsfc.nasa.gov/planetary/factsheet/
 *
 * Calculations are done in real scale and scaled to fit canvas. Need to find a better way of
 * doing that, also implement zoom.
 *
 */
package orbits2;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
//import static java.lang.Math.max;
//import static java.lang.Math.abs;
import java.util.ArrayList;

/**
 *
 * @author raf
 */
public class Orbits2 implements KeyListener, 
                               MouseMotionListener, 
                               MouseWheelListener {

    MyDrawPanel panel;
    ArrayList<Planet> planets = new ArrayList<>();
    
    boolean debug = false;
    boolean paused = false;
    int timer = 10;
    
    boolean showOrbits = false;
    
    double scale;
    int shiftX, shiftY;
    
    public static void main(String[] args) {
        System.out.println("main()");
        Orbits2 orbit = new Orbits2();
        orbit.setup();
    }
    
    public void setup() {
        System.out.println("setup()");
        JFrame frame = new JFrame("Orbits");
        panel = new MyDrawPanel();
        frame.getContentPane().add(panel);
        frame.addKeyListener(this);
        frame.addMouseMotionListener(this);
        frame.addMouseWheelListener(this);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1000, 1000);
        frame.setResizable(true);
        frame.setVisible(true);
        initPlanets();
        start();
    }
    
    private void initPlanets() {
        System.out.println("initPlanets()");
        
        rescale();
        

    }
    private void addPlanet() {
        planets.add(new Planet("nameless",
                               pow(scale, 3), 
                               scale*10,
                               new double[] {0.0, 0.0}, 
                               new double[] {(mouseX-shiftX)*scale,
                                             (mouseY-shiftY)*scale},
                               Color.blue));
    }
    
    private void rescale() {
        int w = panel.getWidth();
        int h = panel.getHeight();
        scale = 1.0;
        shiftX = w/2;
        shiftY = h/2;
    }
    
    private void start() {
        System.out.println("start()");
        System.out.println(panel.getWidth() + " " + panel.getHeight());
        while (true) {
            
            if (!paused) {
                // First update the speed of all planets
                for(Planet p: planets)
                    p.updateVelocity(planets);
                // and then update their positions
                for(Planet p: planets)
                    p.updateOrbit();
            }
            panel.repaint();
            try {
                Thread.sleep(timer);
            } catch (InterruptedException ex) {
                System.out.println("OW NO! WHAT NOW??");
            }
        }
    }
    
    // KeyPressed interface
    @Override
    public void keyPressed(KeyEvent ev) {
        System.out.print("Key pressed: ");
        switch (ev.getKeyCode()) {
            case KeyEvent.VK_A:
                System.out.println("A");
                System.out.println("   add planet");
                addPlanet();
                break;
            case KeyEvent.VK_R:
                System.out.println("R");
                System.out.println("   remove last added planet");
                if (planets.size() > 0)
                    planets.remove(planets.size()-1);
                break;
            case KeyEvent.VK_O:
                System.out.println("O");
                showOrbits = !showOrbits;
                System.out.println("   showOrbit " + showOrbits);
                panel.repaint();
                break;
            case KeyEvent.VK_E:
                System.out.println("E");
                System.out.println("   erase orbits ");
                for (Planet p: planets)
                    p.orbitPoints = 0;
                panel.repaint();
                break;
            case KeyEvent.VK_P:
                System.out.println("P");
                paused = !paused;
                System.out.println("   paused " + paused);
                break;
            default:
                System.out.println("not implemented");
                break;
        }
    }
    @Override public void keyReleased(KeyEvent ev) {}
    @Override public void keyTyped(KeyEvent ev) {}
    
    // MouseMotionListener interface
    int mouseX, mouseY;
    @Override
    public void mouseMoved(MouseEvent e){
        mouseX = e.getPoint().x - 1; 
        mouseY = e.getPoint().y - 24;
    }
    @Override public void mouseDragged(MouseEvent e) {}
   
    // MouseWheelInterface
    @Override public void mouseWheelMoved(MouseWheelEvent e) {
       int notches = e.getWheelRotation();
       if (notches > 0)
           scale *= 1.1;
       else if (notches < 0)
           scale /= 1.1;
       shiftX = (2*shiftX + mouseX) / 3;
       shiftY = (2*shiftY + mouseY) / 3;
       System.out.println("Rescale " + scale);
       System.out.println(" Shift " + shiftX + " " + shiftY);
    }
    
    // drawings
    class MyDrawPanel extends JPanel {
        @Override
        public void paintComponent(Graphics gfx) {
            int w = this.getWidth();
            int h = this.getHeight();
            
            gfx.fillRect(0, 0, w, h);
            gfx.setColor(Color.white);
            gfx.drawString("Mouse " + mouseX + " " + mouseY, 10, 10);

            for (Planet p: planets) {
                if (debug) {
                    System.out.println("Drawing" + p.name);
                    System.out.println("   pos " + p.position[0] + " " + p.position[1]);
                    System.out.println("   rad " + p.radius);
                }
                gfx.setColor(p.color);
                int x = (int)((p.position[0]-p.radius)/scale + shiftX);
                int y = (int)((p.position[1]-p.radius)/scale + shiftY);
                int r = (int)(p.radius/scale);
                gfx.fillOval(x, y, 2*r, 2*r);
            }

            if (showOrbits) {
                //if (debug) System.out.println("Showing orbits");
                gfx.setColor(Color.gray);
                for (Planet p: planets) {
                    //if (debug) System.out.println("  " + p.name);
                    for (int i=0; i < p.orbitPoints-1; i++) {
                        //if (debug)
                        //    System.out.println("   " + p.orbit[i][0] + "   " + p.orbit[i][1]);
                        gfx.drawLine((int)(p.orbit[i  ][0]/scale + shiftX), 
                                     (int)(p.orbit[i  ][1]/scale + shiftY),
                                     (int)(p.orbit[i+1][0]/scale + shiftX),
                                     (int)(p.orbit[i+1][1]/scale + shiftY));
                    }
                }
            }
        }
    }
}


class Planet {
    //static final double G = 6.674e-11;
    static final double G = 6.674e-11;
    static double dt = 24*3600;
    
    boolean debug = false;
    
    String name;
    double mass, radius;
    double[] velocity = new double[2];
    double[] position = new double[2];
    final int size = 5000;
    double[][] orbit;
    int orbitPoints;
    Color color;
    
    double[] u = new double[2];
    double distance;
    
    public Planet(String nm, double m, double r,
                  double[] v, double[] p, Color c) {
        System.out.println("Planet() " + nm);
        name = nm;
        mass = m;
        radius = r;
        position = p;
        velocity = v;
        color = c;
        
        orbit = new double[size][2];
        orbitPoints = 0;
        orbit[orbitPoints][0] = position[0];
        orbit[orbitPoints][1] = position[1];
        orbitPoints++;
        System.out.println("   Position " + position[0] + " "+ position[1]);
        System.out.println("   Mass " + mass);
        System.out.println("   Radius " + radius);
    }
    
    public void updateVelocity(ArrayList<Planet> planets) {
        if (debug)
            System.out.println("Updating velocity of " + name);
        double a;
        for (Planet p: planets) {
            if (p != this) {
                if (debug)
                    System.out.println("   with " + p.name);
                
                distance = sqrt(pow(position[0] - p.position[0], 2) + 
                                   pow(position[1] - p.position[1], 2));
        
                // Unit vector pointing towards the other planet
                u[0] = -(position[0] - p.position[0]) / distance;
                u[1] = -(position[1] - p.position[1]) / distance;
        
                // Classic Newton mechanic:
                // F = m1.a = G.m1.m2.a/r² ==> a = G.m2/r²
                a = G * p.mass / pow(distance, 2);
        
                velocity[0] += a*u[0]*dt;
                velocity[1] += a*u[1]*dt;
            }
        }
    }
    public void updateOrbit() {
        
        position[0] += velocity[0]*dt; 
        position[1] += velocity[1]*dt;
        
        // Check if need to re-allocate orbit arrays
        if (orbitPoints == orbit.length) {
            if (debug) System.out.println("Increasing orbit array size...");
            double[][] tmp = new double[orbit.length + size][2];
            System.arraycopy(orbit, 0, tmp, 0, orbit.length);
            orbit = tmp; 
        }
        
        orbit[orbitPoints][0] = position[0];
        orbit[orbitPoints][1] = position[1];
        orbitPoints++;
        
        if (debug) {
            System.out.println("Updating orbit of " + name);
            System.out.println("   Dist " + distance);
            System.out.println("   u " + u[0] + " " + u[1]);
            System.out.println("   vel " + velocity[0] + " " + velocity[1]);
            System.out.println("   pos " + position[0] + " " + position[1]);
        }
    }
}
