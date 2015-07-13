""" This module provides functionality to extend mpld3 
"""
import mpld3, mpld3.plugins, mpld3.utils

class ClickLink(mpld3.plugins.PluginBase):
    """Plugin for following a link on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clicklink", ClickLink);
    ClickLink.prototype = Object.create(mpld3.Plugin.prototype);
    ClickLink.prototype.constructor = ClickLink;
    ClickLink.prototype.requiredProps = ["id"];
    ClickLink.prototype.defaultProps = {
        links: null
    }
    function ClickLink(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    ClickLink.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        var links = this.props.links;
        
        obj.elements().on("mousedown",
                          function(d, i){ 
                                           window.open(links[i]);
                                        }
                          );
    }
    """
    def __init__(self, points, links):
        self.dict_ = {"type": "clicklink",
                      "id": mpld3.utils.get_id(points),
                      "links": links,
                      }
                 
class MPLSlide(mpld3.plugins.PluginBase):
    JAVASCRIPT = """
         mpld3.Axes.prototype.zoomed = function(propagate) {
            propagate = typeof propagate == "undefined" ? true : propagate;
            if (propagate) {
              var dt0 = this.zoom.translate()[0] - this.zoom.last_t[0];
              var dt1 = this.zoom.translate()[1] - this.zoom.last_t[1];
              var ds = this.zoom.scale() / this.zoom.last_s;
              this.zoom_x.translate([ this.zoom_x.translate()[0] + dt0, 0 ]);
              this.zoom_x.scale(this.zoom_x.scale() * ds);

              this.zoom.last_t = this.zoom.translate();
              this.zoom.last_s = this.zoom.scale();
              this.sharex.forEach(function(ax) {
                ax.zoom_x.translate(this.zoom_x.translate()).scale(this.zoom_x.scale());
              }.bind(this));

              this.sharex.forEach(function(ax) {
                ax.zoomed(false);
              });
            }
            for (var i = 0; i < this.elements.length; i++) {
              this.elements[i].zoomed();
            }
          };
        
            mpld3.ZoomPlugin = mpld3_ZoomPlugin;
            mpld3.register_plugin("zoom", mpld3_ZoomPlugin);
            mpld3_ZoomPlugin.prototype = Object.create(mpld3.Plugin.prototype);
            mpld3_ZoomPlugin.prototype.constructor = mpld3_ZoomPlugin;
            mpld3_ZoomPlugin.prototype.requiredProps = [];
            mpld3_ZoomPlugin.prototype.defaultProps = {
                button: true,
                enabled: null
            };
            function mpld3_ZoomPlugin(fig, props) {
                mpld3.Plugin.call(this, fig, props);
                if (this.props.enabled === null) {
                    this.props.enabled = !this.props.button;
                }
                var enabled = this.props.enabled;
                if (this.props.button) {
                    var ZoomButton = mpld3.ButtonFactory({
                        buttonID: "zoom",
                        sticky: true,
                        actions: [ "scroll", "drag" ],
                        onActivate: this.activate.bind(this),
                        onDeactivate: this.deactivate.bind(this),
                        onDraw: function() {
                            this.setState(enabled);
                        },
                        icon: function() {
                            return mpld3.icons["move"];
                        }
                    });
                this.fig.buttons.push(ZoomButton);
                }
            }
            mpld3_ZoomPlugin.prototype.activate = function() {
                this.fig.enable_zoom();
            };
            mpld3_ZoomPlugin.prototype.deactivate = function() {
                this.fig.disable_zoom();
            };
            mpld3_ZoomPlugin.prototype.draw = function() {
                if (this.props.enabled) this.fig.enable_zoom(); else this.fig.disable_zoom();
            };
        """   
    def __init__(self, button=True, enabled=None):
        if enabled is None:
            enabled = not button
        self.dict_ = {"type": "zoom",
                      "button": button,
                      "enabled": enabled}
                      
class Tooltip(mpld3.plugins.PointHTMLTooltip):
    JAVASCRIPT = ""
    def __init__(self, points, labels=None,
                 hoffset=0, voffset=10, css=None):
        super(Tooltip, self).__init__(points, labels, hoffset, voffset, "")

class LineTooltip(mpld3.plugins.LineHTMLTooltip):
    JAVASCRIPT = ""
    def __init__(self, line, label=None, hoffset=0, voffset=10, css=None):
        super(LineTooltip, self).__init__(line, label, hoffset, voffset, "")
