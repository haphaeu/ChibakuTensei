/*
 * ChibakuTensei
 *
 * Planets orbiting around each other, or colliding =)
 *
 * Calculations are done using classic Newton mechanics.
 *
 *
 */
package chibakutensei;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author raf
 */
public class ChibakuTensei implements KeyListener, 
                                MouseListener,
                                MouseMotionListener, 
                                MouseWheelListener,
                                Runnable {

    MyDrawPanel panel;
    ArrayList<Planet> planets = new ArrayList<>();
    
    boolean debug = false;
    boolean paused = false;
    boolean showTimers = true;
    int timer = 50; // ms
    long updateOrbitTime, updateVelocityTime, checkCollisionsTime, repaintTime; //ns
    int proc_time;  // ms, time required to process and repaint
    int frame_time;  // ms, time to process, draw a frame
    
    boolean showOrbits = true;
    boolean addingPlanetMode = false;
    
    double scale;
    int shiftX, shiftY;
    
    static Thread one;
    
    public static void main(String[] args) {
        System.out.println("main()");
        ChibakuTensei orbit = new ChibakuTensei();

        one = new Thread(orbit);
        one.setName("one");
        orbit.setup();
    }
    
    public void setup() {
        System.out.println("setup()");
        JFrame frame = new JFrame("Orbits");
        panel = new MyDrawPanel();
        frame.getContentPane().add(panel);
        frame.addKeyListener(this);
        frame.addMouseListener(this);
        frame.addMouseMotionListener(this);
        frame.addMouseWheelListener(this);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800, 600);
        frame.setResizable(true);
        frame.setVisible(true);
        rescale();
        initPlanets();
        paused = true;
        start();
    }
    
    private void initPlanets() {
        System.out.println("initPlanets()");
        
        Random rand = new Random();
        int w = panel.getWidth();
        int h = panel.getHeight();
        int numberPlanets = 500;
        double r, m, rho;
        double pad = 150.0;
        double pxmax = 1*w - 2*pad;
        double pymax = 1*h - 2*pad;
        double vmin = -1.0;
        double vmax = 1.0;
        
        for (int i=0; i < numberPlanets; i++) {
            r = rand.nextDouble() * 2.0 + 1.0;
            rho = rand.nextDouble() * 0.1 + 0.01;
            m = 4.0 * rho * pow(r, 3);
            planets.add(new Planet("nameless", 
                    pow(scale, 3) * m,  //  mass
                    scale * r,          // radius
                    new double[] {scale * (rand.nextDouble() * (vmax - vmin) + vmin),
                                  scale * (rand.nextDouble() * (vmax - vmin) + vmin)},
                    new double[] {scale * (pxmax * rand.nextDouble() + pad - shiftX),
                                  scale * (pymax * rand.nextDouble() + pad - shiftY)},
                    new Color(rand.nextFloat(), rand.nextFloat(), rand.nextFloat())));
        }
        //rescale();
    }
    
    private void addPlanet() {
        planets.add(new Planet("nameless",
                               100*pow(scale, 3), 
                               scale*10,
                               new double[] {velX, velY}, 
                               new double[] {(newPlanetX-shiftX)*scale,
                                             (newPlanetY-shiftY)*scale},
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
        long t0, t1, t2, t3;
        int i;
        Planet p;
        
        one.start();
        
        while (true) {
            t0 = System.nanoTime();
            if (!paused) {
                
                // First update the speed of all planets
                for (i=0; i < planets.size(); i++) {
                    p = planets.get(i);
                    p.updateVelocity(planets);
                }
                t1 = System.nanoTime();
                
                // and then update their positions
                
                for (i=0; i < planets.size(); i++) {
                    try {
                        p = planets.get(i);
                        p.updateOrbit();
                    } catch (IndexOutOfBoundsException ex) {
                        ;
                    }
                    
                }
                
                t2 = System.nanoTime();
                checkCollisions();
                t3 = System.nanoTime();
                
                updateVelocityTime = t1 - t0;
                updateOrbitTime = t2 - t1;
                checkCollisionsTime = t3 - t2;
            }
            
            panel.repaint();
            
            proc_time = (int)((System.nanoTime() - t0) / 1e6);
            
            try {
                Thread.sleep(timer - Math.min(timer, proc_time));
            } catch (InterruptedException ex) {
                System.out.println("OW NO! WHAT NOW??");
            }
        }
    }
    private synchronized void checkCollisions() {
        
        one.notify();
        System.out.println(one.getState());

        
    }
    
    @Override
    public synchronized void run() {
        while (true) {
            
            try {
                one.wait();
            } catch (InterruptedException ex) {
                Logger.getLogger(ChibakuTensei.class.getName()).log(Level.SEVERE, null, ex);
            }
            
            for  (int i=0; i < planets.size()-1; i++) {
                Planet pi = planets.get(i);
                for (int j=i+1; j<planets.size(); j++) {
                    Planet pj = planets.get(j);
                    double dist = sqrt(pow(pi.position[0] - pj.position[0], 2) +
                                       pow(pi.position[1] - pj.position[1], 2));
                    if (dist < pi.radius + pj.radius) {
                        //System.out.println("Chibaku tensei!");
                        pi.velocity[0] = (pi.velocity[0]*pi.mass + pj.velocity[0] * pj.mass) / (pi.mass + pj.mass);
                        pi.velocity[1] = (pi.velocity[1]*pi.mass + pj.velocity[1] * pj.mass) / (pi.mass + pj.mass);
                        pi.mass += pj.mass;
                        pi.radius = pow(pow(pi.radius, 3) + pow(pj.radius, 3), 0.33333333);
                        pi.orbitPoints = 0;
                        planets.remove(pj);
                    }
                }
            }
        }
    }
    
    // KeyPressed interface
    int newPlanetX, newPlanetY;
    @Override
    public void keyPressed(KeyEvent ev) {
        System.out.print("Key pressed: ");
        switch (ev.getKeyCode()) {
            case KeyEvent.VK_A:
                if (!addingPlanetMode) {
                    System.out.println("A");
                    System.out.println("   adding planet");
                    newPlanetX = mouseX;
                    newPlanetY = mouseY;
                    newPlanetSpeedX = mouseX;
                    newPlanetSpeedY = mouseY;
                    addingPlanetMode = true;
                    paused = true;
                }
                break;
            case KeyEvent.VK_I:
                System.out.println("I");
                System.out.println("   re-Initi random planets");
                boolean pause_state = paused;
                paused = true;
                initPlanets();
                paused = pause_state;
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
                planets.forEach((p) -> { p.orbitPoints = 0; });
                panel.repaint();
                break;
            case KeyEvent.VK_P:
                System.out.println("P");
                paused = !paused;
                System.out.println("   paused " + paused);
                break;
            case KeyEvent.VK_T:
                System.out.println("T");
                showTimers = !showTimers;
                System.out.println("   showTimers " + showTimers);
                break;
            case KeyEvent.VK_UP:
                System.out.println("up");
                System.out.println("   scroll up ");
                shiftY -= 10;
                break;
            case KeyEvent.VK_DOWN:
                System.out.println("down");
                System.out.println("   scroll down ");
                shiftY += 10;
                break;
            case KeyEvent.VK_LEFT:
                System.out.println("left");
                System.out.println("   scroll left ");
                shiftX -= 10;
                break;
            case KeyEvent.VK_RIGHT:
                System.out.println("right");
                System.out.println("   scroll right ");
                shiftX += 10;
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
    int newPlanetSpeedX, newPlanetSpeedY;
    double velX, velY;
    @Override
    public void mouseMoved(MouseEvent e){
        mouseX = e.getPoint().x - 1; 
        mouseY = e.getPoint().y - 24;
        
        if (addingPlanetMode) {
            newPlanetSpeedX = mouseX;
            newPlanetSpeedY = mouseY;
            if (abs(newPlanetSpeedX - newPlanetX) >= 10 ||
                abs(newPlanetSpeedY - newPlanetY) >= 10) {
                velX = scale * (double)(newPlanetSpeedX - newPlanetX) / 50.0;
                velY = scale * (double)(newPlanetSpeedY - newPlanetY) / 50.0;
            } else {
                velX = 0;
                velY = 0;
            }
        }
    }
   
    // MouseWheelInterface
    @Override public void mouseWheelMoved(MouseWheelEvent e) {
        if (!addingPlanetMode) {
            int cX = (int)((mouseX-shiftX)*scale);
            int cY = (int)((mouseY-shiftY)*scale);
            int notches = e.getWheelRotation();
            if (notches > 0)
                scale *= 1.1;
            else if (notches < 0)
                scale /= 1.1;
            shiftX = mouseX - (int)(cX/scale);
            shiftY = mouseY - (int)(cY/scale);
            System.out.println("Rescale " + scale);
            System.out.println(" Shift " + shiftX + " " + shiftY);
        }
    }

    // MouseListener
    
    @Override public void mouseClicked(MouseEvent me) {
        if (addingPlanetMode){
            addingPlanetMode = false;
            paused = false;
            addPlanet();
        }
    }

    int panX, panY;
    boolean panMode = false;
    @Override public void mousePressed(MouseEvent me) {
        if (!addingPlanetMode) {
            //panMode = true;
            panX = me.getX();
            panY = me.getY();
        }
    }
    @Override public void mouseDragged(MouseEvent e) {
        //if (panMode) {
            shiftX += e.getX() - panX;
            shiftY += e.getY() - panY;
            panX = e.getX();
            panY = e.getY();
        //}
    }
    @Override public void mouseReleased(MouseEvent me) {
        //if (panMode) {
        //    panMode = false;
        //}
    }

    @Override public void mouseEntered(MouseEvent me) {}

    @Override public void mouseExited(MouseEvent me) {}
    
    // drawings
    class MyDrawPanel extends JPanel {
        @Override
        public void paintComponent(Graphics gfx) {
            long t0 = System.nanoTime();
            int w = this.getWidth();
            int h = this.getHeight();
            
            gfx.fillRect(0, 0, w, h);

            Planet p;
            for (int j=0; j < planets.size(); j++) {
                p = planets.get(j);
                        
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
                for (int j=0; j<planets.size(); j++) {
                    p = planets.get(j);
                    //if (debug) System.out.println("  " + p.name);
                    //for (int i=Math.max(0, p.orbitPoints-1000); i < p.orbitPoints-10; i += 10) {
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
            
            gfx.setColor(Color.white);
            gfx.drawString(String.format("Scale %.1f", scale), 10, 10);
            gfx.drawString(String.format("Planets %d", planets.size()), 10, 25);
            if (showTimers) {
                gfx.drawString("Mouse " + mouseX + " " + mouseY, 10, h-10);
                gfx.drawString(String.format("fps %.1f", 1000.0/(proc_time+timer)), 10, h - 25);
                gfx.drawString(String.format("update velocity %5.1fms", updateVelocityTime/1e6), 10, h - 40);
                gfx.drawString(String.format("update orbits %5.1fms", updateOrbitTime/1e6), 10, h - 55);
                gfx.drawString(String.format("check collisions %5.1fms", checkCollisionsTime/1e6), 10, h - 70);
                gfx.drawString(String.format("repaint %5.1fms", repaintTime/1e6), 10, h - 85);
            }
            if (addingPlanetMode) {
                gfx.setColor(Color.gray);
                gfx.drawString(String.format("Velocity %.2f %.2f", velX, velY), 10, 40);
                gfx.drawLine(newPlanetX, newPlanetY,
                             newPlanetSpeedX, newPlanetSpeedY);
            }
            
            repaintTime = System.nanoTime() - t0;
        }
    }
}


class Planet {
    //static final double G = 6.674e-11;
    static final double G = 1.0;
    static double dt = 1;
    
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
        //System.out.println("new Planet() ");
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
        //System.out.println("   Position " + position[0] + " "+ position[1]);
        //System.out.println("   Mass " + mass);
        //System.out.println("   Radius " + radius);
        //System.out.println("   Speed " + velocity[0] + " " + velocity[1]);
    }
    
    public void updateVelocity(ArrayList<Planet> planets) {
        double a;
        Planet p;
        for (int i=0; i < planets.size(); i++) {
            p = planets.get(i);
            if (p != this) {
                distance = sqrt(pow(position[0] - p.position[0], 2) + 
                                   pow(position[1] - p.position[1], 2));
        
                // Unit vector pointing towards the other planet
                u[0] = -(position[0] - p.position[0]) / distance;
                u[1] = -(position[1] - p.position[1]) / distance;
        
                // Classic Newton mechanic:
                // F = m1.a = G.m1.m2/r² ==> a = G.m2/r²
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

