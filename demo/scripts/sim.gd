extends PanelContainer

@onready var result_rect := %ResultRect

var sim : FluidSim = FluidSim.new()
var rel = Vector2.ZERO

func _ready():
	sim.dimension = 256
	sim.diffusion = 0.00001
	sim.init()

func _input(event):
	if event is InputEventMouseMotion:
		rel = event.relative

func _process(delta):
	if Input.is_mouse_button_pressed(MOUSE_BUTTON_LEFT):
		var v = Vector2i(result_rect.get_local_mouse_position() * (float(sim.dimension) / result_rect.size.x))
		for y in range(-20, 21):
			for x in range(-20, 21):
				if Vector2(v).distance_to(Vector2(v) + Vector2(x, y)) < 5.0:
					sim.add_density(v + Vector2i(x, y), 0.25)
		sim.add_velocity(v, rel)

func _physics_process(delta):
	sim.step(delta)
	result_rect.texture = sim.get_density_texture()
