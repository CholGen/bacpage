from textual.app import ComposeResult
from textual.message import Message
from textual.validation import Number
from textual.widgets import Button, Static, Input
from messages import DEFAULT_INFORMATION, MEAN_MESSAGE


class InfoButton( Button ):
    class Description( Message ):

        def __init__( self, description: str ) -> None:
            self.description = description
            super().__init__()

    def __init__( self, label: str, description: str ) -> None:
        self.description = description
        super().__init__( label=label )

    def on_enter( self ) -> None:
        self.post_message( self.Description( self.description ) )

    def on_leave( self ) -> None:
        self.post_message( self.Description( DEFAULT_INFORMATION ) )


class InfoStatic( Static ):
    class Hover( Message ):

        def __init__( self ) -> None:
            super().__init__()

    def on_enter( self ) -> None:
        self.post_message( self.Hover() )


class HoverInput( Input ):
    class Hover( Message ):

        def __init__( self ) -> None:
            super().__init__()

    def on_enter( self ) -> None:
        self.post_message( self.Hover() )


class InfoInput( Static ):
    class Description( Message ):

        def __init__( self, description: str ) -> None:
            self.description = description
            super().__init__()

    def compose( self ) -> ComposeResult:
        yield InfoStatic( "Parameter #1", id="parameterlabel" )
        yield HoverInput( placeholder="2321", validators=[Number( minimum=0, maximum=100 )], id="parameterinput" )

    def __init__( self, description: str ) -> None:
        self.description = description
        super().__init__()

    def on_enter( self ) -> None:
        self.add_class( "hovered" )
        self.post_message( self.Description( self.description ) )

    def on_leave( self ) -> None:
        self.remove_class( "hovered" )
        self.post_message( self.Description( DEFAULT_INFORMATION ) )

    def on_info_static_hover( self ) -> None:
        self.add_class( "hovered" )
        self.post_message( self.Description( self.description ) )

    def on_hover_input_hover( self ) -> None:
        self.add_class( "hovered" )
        self.post_message( self.Description( self.description ) )
